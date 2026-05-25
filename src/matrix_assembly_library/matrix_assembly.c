//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../domains_library/mesh_info_data.h"
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"
#include "../utils/utils.h"

#include "assembly_common.c"

#if defined(COMPILE_CUDA) || defined(COMPILE_SYCL)
#define COMPILE_GPU
#endif

#ifdef COMPILE_GPU
#include "../gpu_utils/accel_utils.h"
#endif

INIT_ASSEMBLY_MATRIX(set_initial_conditions_fvm) {

    real_cpu alpha;

    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    struct cell_node **ac = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;

    OMP(parallel for private(alpha))
    for(uint32_t i = 0; i < active_cells; i++) {
        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v * alpha;
    }
}

INIT_ASSEMBLY_MATRIX(set_initial_conditions_from_odes) {

    real_cpu alpha;

    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    struct cell_node **active_volumes = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;

    int n_equations_cell_model = the_ode_solver->model_data.number_of_ode_equations;
    real *sv = the_ode_solver->sv;

#ifdef COMPILE_GPU
    real *vms = NULL;
    size_t mem_size = the_ode_solver->original_num_cells * sizeof(real);

    if(the_ode_solver->gpu) {
        vms = MALLOC_BYTES(real, mem_size);
        // check_cuda_error(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));
        memcpy_device(vms, sv, mem_size, DEVICE_TO_HOST);
    }
#endif

    OMP(parallel for private(alpha))
    for(uint32_t i = 0; i < active_cells; i++) {
        alpha = ALPHA(beta, cm, dt, active_volumes[i]->discretization.x, active_volumes[i]->discretization.y, active_volumes[i]->discretization.z);

        if(the_ode_solver->gpu) {
#ifdef COMPILE_GPU
            active_volumes[i]->v = vms[active_volumes[i]->sv_position];
            active_volumes[i]->b = vms[active_volumes[i]->sv_position] * alpha;
#endif
        } else {
            active_volumes[i]->v = sv[active_volumes[i]->sv_position * n_equations_cell_model];
            active_volumes[i]->b = sv[active_volumes[i]->sv_position * n_equations_cell_model] * alpha;
        }
    }
#ifdef COMPILE_GPU
    free(vms);
#endif
}

ASSEMBLY_MATRIX(random_sigma_discretization_matrix) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    srand((unsigned int)time(NULL));

    real_cpu modifiers[4] = {0.0f, 0.1f, 0.5f, 1.0f};

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real_cpu r;

        OMP(critical)
        r = modifiers[rand_range(4)];

        real sigma_x_new = sigma_x * r;
        real sigma_y_new = sigma_y * r;
        real sigma_z_new = sigma_z * r;

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

ASSEMBLY_MATRIX(source_sink_discretization_matrix) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    real channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_width, config, "channel_width");

    real channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_length, config, "channel_length");

    bool inside;

    // real side_length_x = the_grid->mesh_side_length.x;
    real side_length_y = the_grid->mesh_side_length.y;
    // real side_length_z = the_grid->mesh_side_length.z;

    real region_height = (side_length_y - channel_width) / 2.0;

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real sigma_x_new;
        real sigma_y_new;
        real sigma_z_new;

        double x = ac[i]->center.x;
        double y = ac[i]->center.y;
        //        double z = ac[i]->center.z;

        // Check region 1
        inside = (x >= 0.0) && (x <= channel_length) && (y >= 0.0) && (y <= region_height);

        // Check region 2
        inside |= (x >= 0.0) && (x <= channel_length) && (y >= region_height + channel_width) && (y <= side_length_y);

        if(inside) {
            sigma_x_new = 0.0;
            sigma_y_new = 0.0;
            sigma_z_new = 0.0;
        } else {
            sigma_x_new = sigma_x;
            sigma_y_new = sigma_y;
            sigma_z_new = sigma_z;
        }

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;
    }

    // Then, we fill the discretization matrix
    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

ASSEMBLY_MATRIX(source_sink_discretization_matrix_with_different_sigma) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    real channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_width, config, "channel_width");

    real channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_length, config, "channel_length");

    real source_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, source_factor, config, "source_factor");

    real sink_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sink_factor, config, "sink_factor");

    bool inside_3, inside_4;

    real side_length_x = the_grid->mesh_side_length.x;
    real side_length_y = the_grid->mesh_side_length.y;
    // real side_length_z = the_grid->mesh_side_length.z;

    real region_height = (side_length_y - channel_width) / 2.0;

    // Set the conductivities for each cell on the grid
    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real sigma_x_new = sigma_x;
        real sigma_y_new = sigma_y;
        real sigma_z_new = sigma_z;

        real x = ac[i]->center.x;
        real y = ac[i]->center.y;
        //        real z = ac[i]->center.z;

        // Check region 3
        inside_3 = (x >= 0.0) && (x < channel_length) && (y >= region_height) && (y <= region_height + channel_width);

        if(inside_3) {
            sigma_x_new = sigma_x * source_factor;
            sigma_y_new = sigma_y * source_factor;
            sigma_z_new = sigma_z * source_factor;
        }

        // Check region 4
        inside_4 = (x >= channel_length) && (x <= side_length_x) && (y >= 0.0) && (y <= side_length_y);

        if(inside_4) {
            sigma_x_new = sigma_x * sink_factor;
            sigma_y_new = sigma_y * sink_factor;
            sigma_z_new = sigma_z * sink_factor;
        }

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;
    }

    // Then, we fill the discretization matrix
    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

ASSEMBLY_MATRIX(homogeneous_sigma_assembly_matrix) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    if(!sigma_initialized) {
        OMP(parallel for)
        for(i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
        }

        // sigma_initialized = true;
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

ASSEMBLY_MATRIX(anisotropic_sigma_assembly_matrix) {

    if(the_grid->adaptive) {
        log_error_and_exit("anisotropic_sigma_assembly_matrix function does not support mesh adaptivity yet!. Aborting!\n");
    }

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    //      D tensor    //
    // | sx    sxy   sxz |
    // | sxy   sy    syz |
    // | sxz   syz   sz  |
    real_cpu D[3][3];
    int i;

    real_cpu sigma_l = 0.0;
    real_cpu sigma_t = 0.0;
    real_cpu sigma_n = 0.0;

    char *fiber_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(fiber_file, config, "fibers_file");

    bool fibers_in_mesh = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibers_in_mesh, config, "fibers_in_mesh");

    struct fiber_coords *fibers = NULL;

    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l, config, "sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t, config, "sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n, config, "sigma_n");

    real_cpu *f = NULL;
    real_cpu *s = NULL;
    real_cpu *n = NULL;

    if(fiber_file) {
        log_info("Loading mesh fibers\n");
        fibers = read_fibers(fiber_file, true);
    } else if(!fibers_in_mesh) {
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(f, config, "f", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(s, config, "s", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(n, config, "n", 3);

        if(!f) {
            f = malloc(sizeof(real_cpu) * 3);
            f[0] = 1.0;
            f[1] = 0.0;
            f[2] = 0.0;
        }

        if(!s) {
            s = malloc(sizeof(real_cpu) * 3);
            s[0] = 0.0;
            s[1] = 1.0;
            s[2] = 0.0;
        }

        if(!n) {
            n = malloc(sizeof(real_cpu) * 3);
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 1.0;
        }
    }

    OMP(parallel for private(D))
    for(i = 0; i < num_active_cells; i++) {

        if(fibers) {
            int fiber_index = ac[i]->original_position_in_file;

            if(fiber_index == -1) {
                log_error_and_exit("fiber_index should not be -1, but it is for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y,
                                   ac[i]->center.z);
            }

            if(sigma_t == sigma_n) {
                calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
            } else {
                calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
            }
            ac[i]->sigma.fibers = fibers[fiber_index];
        } else if(fibers_in_mesh) {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, ac[i]->sigma.fibers.f, sigma_l, sigma_t);
            } else {
                calc_tensor(D, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l, sigma_t, sigma_n);
            }

        } else {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, f, sigma_l, sigma_t);
            } else {
                calc_tensor(D, f, s, n, sigma_l, sigma_t, sigma_n);
            }
        }

        ac[i]->sigma.x = D[0][0];
        ac[i]->sigma.y = D[1][1];
        ac[i]->sigma.z = D[2][2];

        ac[i]->sigma.xy = D[0][1];
        ac[i]->sigma.xz = D[0][2];
        ac[i]->sigma.yz = D[1][2];
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(ac[i]);
    }

    free(f);
    free(s);
    free(n);
}

ASSEMBLY_MATRIX(homogeneous_sigma_with_a_factor_assembly_matrix) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config, "sigma_factor");

    if(!sigma_initialized) {
        OMP(parallel for)
        for(i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x * sigma_factor;
            ac[i]->sigma.y = sigma_y * sigma_factor;
            ac[i]->sigma.z = sigma_z * sigma_factor;
        }

        // sigma_initialized = true;
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

ASSEMBLY_MATRIX(fibrotic_region_with_sigma_factor_assembly_matrix) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config, "sigma_factor");

    real min_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, min_x, config, "region_min_x");

    real max_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, max_x, config, "region_max_x");

    real min_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, min_y, config, "region_min_y");

    real max_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, max_y, config, "region_max_y");

    real min_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, min_z, config, "region_min_z");

    real max_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, max_z, config, "region_max_z");

    bool inside;

    // Set the conductivities for each cell on the grid
    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real sigma_x_new = sigma_x;
        real sigma_y_new = sigma_y;
        real sigma_z_new = sigma_z;

        real x = ac[i]->center.x;
        real y = ac[i]->center.y;
        real z = ac[i]->center.z;

        // Check if inside the region
        inside = (x >= min_x) && (x <= max_x) && (y >= min_y) && (y <= max_y) && (z >= min_z) && (z <= max_z);

        if(inside) {
            sigma_x_new = sigma_x * sigma_factor;
            sigma_y_new = sigma_y * sigma_factor;
            sigma_z_new = sigma_z * sigma_factor;
        }

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;
    }

    // Then, we fill the discretization matrix
    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

ASSEMBLY_MATRIX(heterogenous_sigma_with_factor_assembly_matrix) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config, "seed");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config, "sigma_factor");

    log_info("Reducing conductivity from %.2lf %% of cells\n", phi * 100.0);

    // Initialize the seed for the fibrosis
    srand(seed);

    log_info("Using %u as seed\n", seed);

    if(!sigma_initialized) {

        FOR_EACH_CELL(the_grid) {

            if(cell->active) {
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                if(p < phi) {
                    cell->sigma.x = sigma_x * sigma_factor;
                    cell->sigma.y = sigma_y * sigma_factor;
                    cell->sigma.z = sigma_z * sigma_factor;
                } else {
                    cell->sigma.x = sigma_x;
                    cell->sigma.y = sigma_y;
                    cell->sigma.z = sigma_z;
                }
            }
        }

        // sigma_initialized = true;
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

// This function will read the fibrotic regions and for each cell that is inside the region we will
// reduce its conductivity value based on the 'sigma_factor'.
ASSEMBLY_MATRIX(heterogenous_sigma_with_factor_assembly_matrix_from_file) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    char *fib_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(fib_file, config, "fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config, "size");

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config, "sigma_factor");

    if(!sigma_initialized) {
        OMP(parallel for)
        for(uint32_t i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
        }

        // sigma_initialized = true;
    }

    // Reading the fibrotic regions from the input file
    FILE *file = fopen(fib_file, "r");

    if(!file) {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 7);
        if(scar_mesh[i] == NULL) {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    uint32_t i = 0;
    while(fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3], &scar_mesh[i][4],
                 &scar_mesh[i][5], &scar_mesh[i][6]) != EOF) {
        i++;
    }

    fclose(file);

    uint32_t num_fibrotic_regions = i;

    // Pass through all the cells of the grid and check if its center is inside the current
    // fibrotic region
    OMP(parallel for)
    for(uint32_t j = 0; j < num_fibrotic_regions; j++) {

        struct cell_node *grid_cell = the_grid->first_cell;

        real_cpu b_center_x = scar_mesh[j][0];
        real_cpu b_center_y = scar_mesh[j][1];

        real_cpu b_h_dx = scar_mesh[j][3];
        real_cpu b_h_dy = scar_mesh[j][4];

        bool active = (bool)(scar_mesh[j][6]);

        while(grid_cell != 0) {
            if(grid_cell->active) {
                real_cpu center_x = grid_cell->center.x;
                real_cpu center_y = grid_cell->center.y;
                //               real_cpu half_dx = grid_cell->discretization.x/2.0;
                //                real_cpu half_dy = grid_cell->discretization.y/2.0;

                struct point_3d p;
                struct point_3d q;

                p.x = b_center_y + b_h_dy;
                p.y = b_center_y - b_h_dy;

                q.x = b_center_x + b_h_dx;
                q.y = b_center_x - b_h_dx;

                // Check if the current cell is inside the fibrotic region
                if(center_x > q.y && center_x < q.x && center_y > p.y && center_y < p.x) {
                    if(active == 0) {
                        grid_cell->sigma.x = sigma_x * sigma_factor;
                        grid_cell->sigma.y = sigma_y * sigma_factor;
                        grid_cell->sigma.z = sigma_z * sigma_factor;
                    }
                }
            }

            grid_cell = grid_cell->next;
        }
    }

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }

    for(int k = 0; k < fib_size; k++) {
        free(scar_mesh[k]);
    }

    free(scar_mesh);
}

// This function will generate the fibrotic region file for the 120um x 120um grid by reescaling
// the original Scientific Reports 4b grid from 40000um side_length to 48000um
// rescale_factor = 1.2
ASSEMBLY_MATRIX(heterogenous_fibrotic_region_file_write_with_input_file) {

    initialize_diagonal_elements(the_solver, the_grid);

    char *fib_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(fib_file, config, "fibrosis_file");

    char *new_fib_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(new_fib_file, config, "rescaled_fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config, "size");

    real rescale_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, rescale_factor, config, "rescale_factor");

    FILE *file = fopen(fib_file, "r");

    if(!file) {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    // Read and store the original positions of the fibrotic regions
    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 7);
        if(scar_mesh[i] == NULL) {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    for(int i = 0; i < fib_size; i++) {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3], &scar_mesh[i][4],
               &scar_mesh[i][5], &scar_mesh[i][6]);
    }

    fclose(file);

    // Write the new fibrotic region file based on the 'rescale_factor'
    FILE *fileW = fopen(new_fib_file, "w");

    if(!file) {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    // Multiple the positions of each scar by a rescale factor
    for(int i = 0; i < fib_size; i++) {
        scar_mesh[i][0] = scar_mesh[i][0] * rescale_factor;
        scar_mesh[i][1] = scar_mesh[i][1] * rescale_factor;
        scar_mesh[i][2] = scar_mesh[i][2] * rescale_factor;
        scar_mesh[i][3] = scar_mesh[i][3] * rescale_factor;
        scar_mesh[i][4] = scar_mesh[i][4] * rescale_factor;
        scar_mesh[i][5] = scar_mesh[i][5] * rescale_factor;

        fprintf(fileW, "%g,%g,%g,%g,%g,%g,%g\n", scar_mesh[i][0], scar_mesh[i][1], scar_mesh[i][2], scar_mesh[i][3], scar_mesh[i][4], scar_mesh[i][5],
                scar_mesh[i][6]);
    }

    fclose(fileW);

    for(int k = 0; k < fib_size; k++) {
        free(scar_mesh[k]);
    }

    free(scar_mesh);

    // We just leave the program after this ...
    log_info("[!] Finish writing new fibrotic region file '%s'!\n", new_fib_file);
    exit(EXIT_SUCCESS);
}

// This function will generate the fibrotic region file for the 120um x 120um grid by reescaling
// the original Scientific Reports 4b grid from 40000um side_length to 48000um. The fibrotic region
// will be mapped using the same idea used on the 'domains_library' function by using a fixed seed
// for the random number generator.
// rescale_factor = 1.2
ASSEMBLY_MATRIX(heterogenous_fibrotic_region_file_write_using_seed) {

    initialize_diagonal_elements(the_solver, the_grid);

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config, "seed");

    char *new_fib_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(new_fib_file, config, "rescaled_fibrosis_file");

    double rescale_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, rescale_factor, config, "rescale_factor");

    // Write the new fibrotic region file
    FILE *fileW = fopen(new_fib_file, "w+");

    // Initialize the random the generator with the same seed used by the original model
    srand(seed);
    FOR_EACH_CELL(the_grid) {

        if(cell->active) {
            real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
            if(p < phi) {
                // We reescale the cell position using the 'rescale_factor'
                double center_x = cell->center.x * rescale_factor;
                double center_y = cell->center.y * rescale_factor;
                double center_z = cell->center.z * rescale_factor;
                double dx = cell->discretization.x * rescale_factor;
                double dy = cell->discretization.y * rescale_factor;
                double dz = cell->discretization.z * rescale_factor;

                // Then, we write only the fibrotic regions to the output file
                fprintf(fileW, "%g,%g,%g,%g,%g,%g,0\n", center_x, center_y, center_z, dx / 2.0, dy / 2.0, dz / 2.0);
            }
        }
    }

    fclose(fileW);

    // We just leave the program after this ...
    log_info("[!] Finish writing fibrotic region file '%s'!\n", new_fib_file);
    exit(EXIT_SUCCESS);
}

ASSEMBLY_MATRIX(homogeneous_sigma_assembly_matrix_with_fast_endocardium_layer) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    real_cpu fast_endo_layer_scale = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, fast_endo_layer_scale, config, "fast_endo_layer_scale");

    if(!sigma_initialized) {
        OMP(parallel for)
        for(i = 0; i < num_active_cells; i++) {
            // Check if the current cell is tagged as FASTENDO
            real_cpu tag = TISSUE_TYPE(ac[i]);
            if(tag == 0) {
                ac[i]->sigma.x = sigma_x * fast_endo_layer_scale;
                ac[i]->sigma.y = sigma_y * fast_endo_layer_scale;
                ac[i]->sigma.z = sigma_z * fast_endo_layer_scale;
            }
            // Normal type of cell
            else {
                ac[i]->sigma.x = sigma_x;
                ac[i]->sigma.y = sigma_y;
                ac[i]->sigma.z = sigma_z;
            }
        }
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

ASSEMBLY_MATRIX(anisotropic_sigma_assembly_matrix_with_fast_endocardium_layer) {

    if(the_grid->adaptive) {
        log_error_and_exit("anisotropic_sigma_assembly_matrix_with_fast_endocardium_layer function does not support mesh adaptivity yet!. Aborting!\n");
    }

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    //      D tensor    //
    // | sx    sxy   sxz |
    // | sxy   sy    syz |
    // | sxz   syz   sz  |
    real_cpu D[3][3];
    int i;

    real_cpu sigma_l = 0.0;
    real_cpu sigma_t = 0.0;
    real_cpu sigma_n = 0.0;

    char *fiber_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(fiber_file, config, "fibers_file");

    bool fibers_in_mesh = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibers_in_mesh, config, "fibers_in_mesh");

    real_cpu fast_endo_layer_scale = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, fast_endo_layer_scale, config, "fast_endo_layer_scale");

    struct fiber_coords *fibers = NULL;

    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l, config, "sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t, config, "sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n, config, "sigma_n");

    real_cpu *f = NULL;
    real_cpu *s = NULL;
    real_cpu *n = NULL;

    if(fiber_file) {
        log_info("Loading mesh fibers\n");
        fibers = read_fibers(fiber_file, false);
    } else if(!fibers_in_mesh) {
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(f, config, "f", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(s, config, "s", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(n, config, "n", 3);

        if(!f) {
            f = malloc(sizeof(real_cpu) * 3);
            f[0] = 1.0;
            f[1] = 0.0;
            f[2] = 0.0;
        }

        if(!s) {
            s = malloc(sizeof(real_cpu) * 3);
            s[0] = 0.0;
            s[1] = 1.0;
            s[2] = 0.0;
        }

        if(!n) {
            n = malloc(sizeof(real_cpu) * 3);
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 1.0;
        }
    }

    OMP(parallel for private(D))
    for(i = 0; i < num_active_cells; i++) {

        if(fibers) {
            int fiber_index = ac[i]->original_position_in_file;

            if(fiber_index == -1) {
                log_error_and_exit("fiber_index should not be -1, but it is for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y,
                                   ac[i]->center.z);
            }

            if(sigma_t == sigma_n) {
                calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
            } else {
                calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
            }
            ac[i]->sigma.fibers = fibers[fiber_index];
        } else if(fibers_in_mesh) {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, ac[i]->sigma.fibers.f, sigma_l, sigma_t);
            } else {
                calc_tensor(D, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l, sigma_t, sigma_n);
            }

        } else {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, f, sigma_l, sigma_t);
            } else {
                calc_tensor(D, f, s, n, sigma_l, sigma_t, sigma_n);
            }
        }

        // Check if the current cell is tagged as FASTENDO
        // TODO: Try to load the "extra_data" array and check if a cell is fast_endo or not
        real_cpu tag = DTI_MESH_TRANSMURALITY_LABELS(ac[i]);

        if(tag == 3) {
            ac[i]->sigma.x = D[0][0] * fast_endo_layer_scale;
            ac[i]->sigma.y = D[1][1] * fast_endo_layer_scale;
            ac[i]->sigma.z = D[2][2] * fast_endo_layer_scale;

            ac[i]->sigma.xy = D[0][1] * fast_endo_layer_scale;
            ac[i]->sigma.xz = D[0][2] * fast_endo_layer_scale;
            ac[i]->sigma.yz = D[1][2] * fast_endo_layer_scale;
        }
        // Normal type of cell
        else {
            ac[i]->sigma.x = D[0][0];
            ac[i]->sigma.y = D[1][1];
            ac[i]->sigma.z = D[2][2];

            ac[i]->sigma.xy = D[0][1];
            ac[i]->sigma.xz = D[0][2];
            ac[i]->sigma.yz = D[1][2];
        }
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(ac[i]);
    }

    free(f);
    free(s);
    free(n);
}
