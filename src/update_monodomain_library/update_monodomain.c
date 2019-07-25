//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "../alg/grid/grid.h"
#include "../config/update_monodomain_config.h"
#include "../libraries_common/config_helpers.h"
#include "../monodomain/constants.h"
#include "../utils/utils.h"
#include "../single_file_libraries/stb_ds.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif


UPDATE_MONODOMAIN(update_monodomain_default) {

    real_cpu alpha;
    bool use_gpu = the_ode_solver->gpu;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **active_cells = the_grid->active_cells;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt_pde = the_solver->dt;

    int n_equations_cell_model = the_ode_solver->model_data.number_of_ode_equations;
    real *sv = the_ode_solver->sv;

    #ifdef COMPILE_CUDA
    real *vms = NULL;
    size_t mem_size = initial_number_of_cells * sizeof(real);

    if(use_gpu) {
        vms = (real *)malloc(mem_size);
        check_cuda_errors(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));
    }
    #endif

    int i;
    #pragma omp parallel for private(alpha)
    for(i = 0; i < num_active_cells; i++) {
        alpha = ALPHA(beta, cm, dt_pde, active_cells[i]->dx, active_cells[i]->dy, active_cells[i]->dz);

        if(use_gpu) {
            #ifdef COMPILE_CUDA
            active_cells[i]->b = vms[active_cells[i]->sv_position] * alpha;
            #endif
        } else {
            active_cells[i]->b = sv[active_cells[i]->sv_position * n_equations_cell_model] * alpha;
        }
    }
    #ifdef COMPILE_CUDA
    free(vms);
    #endif
}

UPDATE_MONODOMAIN(update_monodomain_ddm) 
{

    real_cpu alpha;
    bool use_gpu = the_ode_solver->gpu;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **active_cells = the_grid->active_cells;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt_pde = the_solver->dt;

    int n_equations_cell_model = the_ode_solver->model_data.number_of_ode_equations;
    real *sv = the_ode_solver->sv;

    //real_cpu kappa_x = the_solver->kappa_x;
    //real_cpu kappa_y = the_solver->kappa_y;
    //real_cpu kappa_z = the_solver->kappa_z;

#ifdef COMPILE_CUDA
    real *vms = NULL;
    size_t mem_size = initial_number_of_cells * sizeof(real);

    if(use_gpu)
    {
        vms = (real *)malloc(mem_size);
        check_cuda_errors(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));
    }
#endif

    int i;

#pragma omp parallel for private(alpha)
    for(i = 0; i < num_active_cells; i++)
    {
        // 1) Calculate alpha for the diagonal element
        alpha = ALPHA(beta, cm, dt_pde, active_cells[i]->dx, active_cells[i]->dy, active_cells[i]->dz);

        if(use_gpu)
        {
#ifdef COMPILE_CUDA
            active_cells[i]->b = vms[active_cells[i]->sv_position] * alpha;
#endif
        }
        else
        {
            active_cells[i]->b = sv[active_cells[i]->sv_position * n_equations_cell_model] * alpha;
        }

        // 2) Calculate kappas
        // We need to capture the neighbours from the current volume
        struct element *cell_elements = active_cells[i]->elements;
        uint32_t max_elements = arrlen(cell_elements);

        real_cpu dx = active_cells[i]->dx;
        real_cpu dy = active_cells[i]->dy;
        real_cpu dz = active_cells[i]->dz;

        real_cpu kappa_x = active_cells[i]->kappa_x;
        real_cpu kappa_y = active_cells[i]->kappa_y;
        real_cpu kappa_z = active_cells[i]->kappa_z;

        for (int j = 1; j < max_elements; j++)
        {
            int k = cell_elements[j].column;

            if (cell_elements[j].direction == 'n') // North cell
            {
                real_cpu multiplier = (dx * dy) / dz;
                if(use_gpu)
                {
#ifdef COMPILE_CUDA
                    active_cells[i]->b -= vms[active_cells[k]->sv_position] * multiplier * kappa_z / dt_pde;
                    active_cells[i]->b += vms[active_cells[i]->sv_position] * multiplier * kappa_z / dt_pde;
#endif
                }
                else
                {
                    active_cells[i]->b -= sv[active_cells[k]->sv_position * n_equations_cell_model] * multiplier * kappa_z / dt_pde;
                    active_cells[i]->b += sv[active_cells[i]->sv_position * n_equations_cell_model] * multiplier * kappa_z / dt_pde;
                }
            }
            else if (cell_elements[j].direction == 's') // South cell
            {
                real_cpu multiplier = (dx * dy) / dz;
                if(use_gpu)
                {
#ifdef COMPILE_CUDA
                    active_cells[i]->b -= vms[active_cells[k]->sv_position] * multiplier * kappa_z / dt_pde;
                    active_cells[i]->b += vms[active_cells[i]->sv_position] * multiplier * kappa_z / dt_pde;
#endif
                }
                else
                {
                    active_cells[i]->b -= sv[active_cells[k]->sv_position * n_equations_cell_model] * multiplier * kappa_z / dt_pde;
                    active_cells[i]->b += sv[active_cells[i]->sv_position * n_equations_cell_model] * multiplier * kappa_z / dt_pde;
                }
            }
            else if (cell_elements[j].direction == 'e') // East cell
            {
                real_cpu multiplier = (dx * dz) / dy;
                if(use_gpu)
                {
#ifdef COMPILE_CUDA
                    active_cells[i]->b -= vms[active_cells[k]->sv_position] * multiplier * kappa_y / dt_pde;
                    active_cells[i]->b += vms[active_cells[i]->sv_position] * multiplier * kappa_y / dt_pde;
#endif
                }
                else
                {
                    active_cells[i]->b -= sv[active_cells[k]->sv_position * n_equations_cell_model] * multiplier * kappa_y / dt_pde;
                    active_cells[i]->b += sv[active_cells[i]->sv_position * n_equations_cell_model] * multiplier * kappa_y / dt_pde;
                }
            }
            else if (cell_elements[j].direction == 'w') // West cell
            {
                real_cpu multiplier = (dx * dz) / dy;
                if(use_gpu)
                {
                    #ifdef COMPILE_CUDA
                    active_cells[i]->b -= vms[active_cells[k]->sv_position] * multiplier * kappa_y / dt_pde;
                    active_cells[i]->b += vms[active_cells[i]->sv_position] * multiplier * kappa_y / dt_pde;
                    #endif
                }
                else
                {
                    active_cells[i]->b -= sv[active_cells[k]->sv_position * n_equations_cell_model] * multiplier * kappa_y / dt_pde;
                    active_cells[i]->b += sv[active_cells[i]->sv_position * n_equations_cell_model] * multiplier * kappa_y / dt_pde;
                }
            }
            else if (cell_elements[j].direction == 'f') // Forward cell
            {
                real_cpu multiplier = (dy * dz) / dx;
                if(use_gpu)
                {
                    #ifdef COMPILE_CUDA
                    active_cells[i]->b -= vms[active_cells[k]->sv_position] * multiplier * kappa_x / dt_pde;
                    active_cells[i]->b += vms[active_cells[i]->sv_position] * multiplier * kappa_x / dt_pde;
                    #endif
                }
                else
                {
                    active_cells[i]->b -= sv[active_cells[k]->sv_position * n_equations_cell_model] * multiplier * kappa_x / dt_pde;
                    active_cells[i]->b += sv[active_cells[i]->sv_position * n_equations_cell_model] * multiplier * kappa_x / dt_pde;
                }
            }
            else if (cell_elements[j].direction == 'b') // Backward cell
            {
                real_cpu multiplier = (dy * dz) / dx;
                if(use_gpu)
                {
                    #ifdef COMPILE_CUDA
                    active_cells[i]->b -= vms[active_cells[k]->sv_position] * multiplier * kappa_x / dt_pde;
                    active_cells[i]->b += vms[active_cells[i]->sv_position] * multiplier * kappa_x / dt_pde;
                    #endif
                }
                else
                {
                    active_cells[i]->b -= sv[active_cells[k]->sv_position * n_equations_cell_model] * multiplier * kappa_x / dt_pde;
                    active_cells[i]->b += sv[active_cells[i]->sv_position * n_equations_cell_model] * multiplier * kappa_x / dt_pde;
                }
            }
        }
    }
#ifdef COMPILE_CUDA
    free(vms);
#endif

}