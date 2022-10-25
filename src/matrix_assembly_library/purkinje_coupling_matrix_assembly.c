//
// Created by bergolho on 29/09/20.
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
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"
#include "../utils/utils.h"
#include "../domains_library/mesh_info_data.h"

#include "assembly_common.c"

INIT_ASSEMBLY_MATRIX(set_initial_conditions_coupling_fvm) {

    real_cpu alpha;

    // Tissue parameters
    struct cell_node **ac = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;

    // Purkinje parameters
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;
    uint32_t active_purkinje_cells = the_grid->purkinje->num_active_purkinje_cells;

    // Common parameters
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;
    uint32_t i;

    // Tissue section
    OMP(parallel for private(alpha))
    for(i = 0; i < active_cells; i++)
    {

        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v * alpha;
    }

    // Purkinje section
    OMP(parallel for private(alpha))
    for(i = 0; i < active_purkinje_cells; i++) {

        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac_purkinje[i]->v = purkinje_initial_v;
        ac_purkinje[i]->b = purkinje_initial_v * alpha;
    }
}

void initialize_diagonal_elements_purkinje (struct monodomain_solver *the_solver, struct grid *the_grid) {

    real_cpu alpha;
    real_cpu dx, dy, dz;

    uint32_t num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac = the_grid->purkinje->purkinje_cells;

    struct node *n = the_grid->purkinje->network->list_nodes;

    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    uint32_t i;

    for (i = 0; i < num_active_cells; i++) {

        dx = ac[i]->discretization.x;
        dy = ac[i]->discretization.y;
        dz = ac[i]->discretization.z;

        alpha = ALPHA(beta, cm, dt, dx, dy, dz);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if (ac[i]->elements != NULL)
            arrfree(ac[i]->elements);

        ac[i]->elements = NULL;
        arrsetcap(ac[i]->elements,n->num_edges);
        arrput(ac[i]->elements, element);

        n = n->next;
    }
}

// For the Purkinje fibers we only need to solve the 1D Monodomain equation
static void fill_discretization_matrix_elements_purkinje (bool has_point_data, real_cpu sigma_x, struct cell_node **grid_cells, uint32_t num_active_cells,
                                                        struct node *pk_node) {

    struct edge *e;
    struct element **cell_elements;

    real_cpu dx, dy, dz;
    real_cpu multiplier;

    int i;

    for (i = 0; i < num_active_cells; i++, pk_node = pk_node->next) {

        cell_elements = &grid_cells[i]->elements;
        dx = grid_cells[i]->discretization.x;
        dy = grid_cells[i]->discretization.x;
        dz = grid_cells[i]->discretization.x;

        multiplier = ((dy * dz) / dx);

        e = pk_node->list_edges;

        // Do the mapping of the edges from the graph to the sparse matrix data structure ...
        while (e != NULL) {

            struct element new_element;

            // Calculate the conductivity between the two neighboring cells
            if (has_point_data) {
                real_cpu sigma_x1 = pk_node->sigma;
                real_cpu sigma_x2 = e->dest->sigma;
                
                if(sigma_x1 != 0.0 && sigma_x2 != 0.0) 
                    sigma_x = (2.0f * sigma_x1 * sigma_x2) / (sigma_x1 + sigma_x2);
            }

            // Neighbour elements ...
            new_element.column = e->id;
            new_element.value = (-sigma_x * multiplier);
            new_element.cell = grid_cells[e->id];

            // Diagonal element ...
            cell_elements[0]->value += (sigma_x * multiplier);

            arrput(grid_cells[i]->elements,new_element);

            e = e->next;
        }
    }
}

ASSEMBLY_MATRIX (purkinje_coupling_assembly_matrix) {

// [TISSUE]
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

    real sigma_purkinje = sigma_x;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real,sigma_purkinje,config,"sigma_purkinje");

    if(!sigma_initialized) {

        OMP(parallel for)
        for (i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
        }

        sigma_initialized = true;
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

// [PURKINJE]
    static bool sigma_purkinje_initialized = false;

    uint32_t num_purkinje_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;
    struct node *pk_node = the_grid->purkinje->network->list_nodes;
    bool has_point_data = the_grid->purkinje->network->has_point_data;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    if(!sigma_purkinje_initialized) {
        // Check if the Purkinje network file has the POINT_DATA section
        if (has_point_data) {
            struct node *tmp = the_grid->purkinje->network->list_nodes;
            uint32_t i = 0;
            while (tmp != NULL)
            {
                // Copy the prescribed conductivity from the Purkinje network file into the ALG Purkinje cell structure
                ac_purkinje[i]->sigma.x = tmp->sigma;

                tmp = tmp->next; i++;
            }
        } 
        // Otherwise, initilize the conductivity of all cells homogenously with the value from the configuration file
        else {
            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++) {
                ac[i]->sigma.x = sigma_purkinje;
            }
        }
        sigma_purkinje_initialized = true;
    }
    fill_discretization_matrix_elements_purkinje(has_point_data,sigma_purkinje,ac_purkinje,num_purkinje_active_cells,pk_node);
}

ASSEMBLY_MATRIX (purkinje_coupling_with_anisotropic_sigma_assembly_matrix) {

// [TISSUE]
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
    real_cpu sigma_purkinje = 0.0;

    char *fiber_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(fiber_file, config, "fibers_file");

    bool fibers_in_mesh = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibers_in_mesh, config, "fibers_in_mesh");


    struct fiber_coords *fibers = NULL;

    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l, config, "sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t, config, "sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n, config, "sigma_n");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real,sigma_purkinje,config,"sigma_purkinje");

    real_cpu *f = NULL;
    real_cpu *s = NULL;
    real_cpu *n = NULL;

    if(fiber_file) {
        log_info("Loading mesh fibers\n");
        fibers = read_fibers(fiber_file, false);
    }
    else if(!fibers_in_mesh) {
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(f, config, "f", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(s, config, "s", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(n, config, "n", 3);

        if(!f) {
            f = malloc(sizeof(real_cpu)*3);
            f[0] = 1.0;
            f[1] = 0.0;
            f[2] = 0.0;
        }

        if(!s) {
            s = malloc(sizeof(real_cpu)*3);
            s[0] = 0.0;
            s[1] = 1.0;
            s[2] = 0.0;
        }

        if(!n) {
            n = malloc(sizeof(real_cpu)*3);
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
                log_error_and_exit("fiber_index should not be -1, but it is for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);
            }

            if(sigma_t == sigma_n) {
                calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
            }
            ac[i]->sigma.fibers = fibers[fiber_index];
        }
        else if(fibers_in_mesh) {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, ac[i]->sigma.fibers.f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l, sigma_t, sigma_n);
            }

        }
        else {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, f, sigma_l, sigma_t);
            }
            else {
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

// [PURKINJE]
    static bool sigma_purkinje_initialized = false;

    uint32_t num_purkinje_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;
    struct node *pk_node = the_grid->purkinje->network->list_nodes;
    bool has_point_data = the_grid->purkinje->network->has_point_data;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    if(!sigma_purkinje_initialized) {
        // Check if the Purkinje network file has the POINT_DATA section
        if (has_point_data) {
            struct node *tmp = the_grid->purkinje->network->list_nodes;
            uint32_t i = 0;
            while (tmp != NULL)
            {
                // Copy the prescribed conductivity from the Purkinje network file into the ALG Purkinje cell structure
                ac_purkinje[i]->sigma.x = tmp->sigma;

                tmp = tmp->next; i++;
            }
        } 
        // Otherwise, initilize the conductivity of all cells homogenously with the value from the configuration file
        else {
            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++) {
                ac[i]->sigma.x = sigma_purkinje;
            }
        }
        sigma_purkinje_initialized = true;
    }
    fill_discretization_matrix_elements_purkinje(has_point_data,sigma_purkinje,ac_purkinje,num_purkinje_active_cells,pk_node);
}

ASSEMBLY_MATRIX (purkinje_coupling_with_anisotropic_sigma_and_fast_endocardium_layer_assembly_matrix) {

// [TISSUE]
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
    real_cpu sigma_purkinje = 0.0;

    char *fiber_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(fiber_file, config, "fibers_file");

    bool fibers_in_mesh = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibers_in_mesh, config, "fibers_in_mesh");

    real_cpu fast_endo_layer_scale = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu ,fast_endo_layer_scale, config, "fast_endo_layer_scale");

    struct fiber_coords *fibers = NULL;

    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l, config, "sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t, config, "sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n, config, "sigma_n");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real,sigma_purkinje,config,"sigma_purkinje");

    real_cpu *f = NULL;
    real_cpu *s = NULL;
    real_cpu *n = NULL;

    if(fiber_file) {
        log_info("Loading mesh fibers\n");
        fibers = read_fibers(fiber_file, false);
    }
    else if(!fibers_in_mesh) {
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(f, config, "f", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(s, config, "s", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(n, config, "n", 3);

        if(!f) {
            f = malloc(sizeof(real_cpu)*3);
            f[0] = 1.0;
            f[1] = 0.0;
            f[2] = 0.0;
        }

        if(!s) {
            s = malloc(sizeof(real_cpu)*3);
            s[0] = 0.0;
            s[1] = 1.0;
            s[2] = 0.0;
        }

        if(!n) {
            n = malloc(sizeof(real_cpu)*3);
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
                log_error_and_exit("fiber_index should not be -1, but it is for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);
            }

            if(sigma_t == sigma_n) {
                calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
            }
            ac[i]->sigma.fibers = fibers[fiber_index];
        }
        else if(fibers_in_mesh) {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, ac[i]->sigma.fibers.f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l, sigma_t, sigma_n);
            }

        }
        else {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, f, s, n, sigma_l, sigma_t, sigma_n);
            }
        }

        // Check if the current cell is tagged as FASTENDO
        real_cpu tag = DTI_MESH_TRANSMURALITY_LABELS(ac[i]);
        if (tag == 3) {
            ac[i]->sigma.x = D[0][0]*fast_endo_layer_scale;
            ac[i]->sigma.y = D[1][1]*fast_endo_layer_scale;
            ac[i]->sigma.z = D[2][2]*fast_endo_layer_scale;

            ac[i]->sigma.xy = D[0][1]*fast_endo_layer_scale;
            ac[i]->sigma.xz = D[0][2]*fast_endo_layer_scale;
            ac[i]->sigma.yz = D[1][2]*fast_endo_layer_scale;
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

// [PURKINJE]
    static bool sigma_purkinje_initialized = false;

    uint32_t num_purkinje_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;
    struct node *pk_node = the_grid->purkinje->network->list_nodes;
    bool has_point_data = the_grid->purkinje->network->has_point_data;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    if(!sigma_purkinje_initialized) {
        // Check if the Purkinje network file has the POINT_DATA section
        if (has_point_data) {
            struct node *tmp = the_grid->purkinje->network->list_nodes;
            uint32_t i = 0;
            while (tmp != NULL)
            {
                // Copy the prescribed conductivity from the Purkinje network file into the ALG Purkinje cell structure
                ac_purkinje[i]->sigma.x = tmp->sigma;

                tmp = tmp->next; i++;
            }
        } 
        // Otherwise, initilize the conductivity of all cells homogenously with the value from the configuration file
        else {
            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++) {
                ac[i]->sigma.x = sigma_purkinje;
            }
        }
        sigma_purkinje_initialized = true;
    }
    fill_discretization_matrix_elements_purkinje(has_point_data,sigma_purkinje,ac_purkinje,num_purkinje_active_cells,pk_node);
}
