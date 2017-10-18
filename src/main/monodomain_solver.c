//
// Created by sachetto on 03/10/17.
//


#include <inttypes.h>
#include <omp.h>
#include <assert.h>
#include <string.h>

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

#include "monodomain_solver.h"
#include "../utils/stop_watch.h"
#include "linear_system_solver.h"
#include "../utils/logfile_utils.h"


static inline double ALPHA (double beta, double cm, double dt, double h) {
    return (((beta * cm) / dt) * UM2_TO_CM2) * pow (h, 3.0);
}

struct monodomain_solver *new_monodomain_solver() {

    struct monodomain_solver *result = (struct monodomain_solver *)malloc (sizeof (struct monodomain_solver));

    result->beta = 0.14;
    result->cm = 1.0;

    return result;
}

void solve_monodomain(struct grid *the_grid, struct monodomain_solver *the_monodomain_solver,
                      struct ode_solver *the_ode_solver,  struct user_options *configs) {

    assert(configs);

    assert(the_grid);
    assert(the_monodomain_solver);
    assert(the_ode_solver);

    struct stim_config_hash *stimuli_configs = configs->stim_configs;
    struct extra_data_config *extra_data_config = configs->extra_data_config;
    struct domain_config * domain_config = configs->domain_config;

    bool has_extra_data = (extra_data_config != NULL);

    int refine_each = the_monodomain_solver->refine_each;
    int derefine_each = the_monodomain_solver->derefine_each;

    double cg_tol = the_monodomain_solver->tolerance;

    int np = the_monodomain_solver->num_threads;

    if (np == 0)
        np = 1;

#if defined(_OPENMP)
    omp_set_num_threads (np);
#endif

    bool redo_matrix;

    bool activity;

    int max_its;

    bool gpu = the_ode_solver->gpu;
    bool jacobi = the_monodomain_solver->use_jacobi;

    int count = 0;

    double refinement_bound = the_monodomain_solver->refinement_bound;
    double derefinement_bound = the_monodomain_solver->derefinement_bound;

    double start_h = domain_config->start_h;
    double max_h = domain_config->max_h;

    bool adaptive = the_grid->adaptive;
    double start_adpt_at = the_monodomain_solver->start_adapting_at;
    bool save_to_file = (configs->out_dir_name != NULL);

    double dt_edp = the_monodomain_solver->dt;
    double finalT = the_monodomain_solver->final_time;

    double beta = the_monodomain_solver->beta;
    double cm = the_monodomain_solver->cm;

    double dt_edo = the_ode_solver->min_dt;

#ifdef COMPILE_CUDA
    if (gpu) {
        int device_count;
        int device = the_ode_solver->gpu_id;
        check_cuda_errors(cudaGetDeviceCount (&device_count));
        struct cudaDeviceProp prop;
        check_cuda_errors(cudaGetDeviceProperties (&prop, the_ode_solver->gpu_id));
        print_to_stdout_and_file ("%d devices available, running on Device %d: %s\n", device_count, device, prop.name);
        check_cuda_errors(cudaSetDevice (device));
    }
#endif

    order_grid_cells (the_grid);
    uint32_t original_num_cells = the_grid->num_active_cells;

    save_old_cell_positions (the_grid);
    update_cells_to_solve (the_grid, the_ode_solver);

    print_to_stdout_and_file ("Setting ODE's initial conditions\n");
    set_ode_initial_conditions_for_all_volumes (the_ode_solver, the_grid->num_active_cells);

    double initial_v = the_ode_solver->model_data.initial_v;

    if (the_monodomain_solver->max_iterations > 0) {
        max_its = the_monodomain_solver->max_iterations;
    } else {
        max_its = (int)the_grid->number_of_cells;
    }

    print_solver_info(the_monodomain_solver, the_ode_solver, the_grid, configs);

    int ode_step = 1;

    if (dt_edp >= dt_edo) {
        ode_step = (int)(dt_edp / dt_edo);
        print_to_stdout_and_file ("Solving EDO %d times before solving PDE\n", ode_step);
    } else {
        print_to_stdout_and_file ("WARNING: EDO time step is greater than PDE time step. Adjusting to EDO time "
                "step: %lf\n",
                dt_edo);
        dt_edp = dt_edo;
    }

    fflush (stdout);

    long ode_total_time = 0, cg_total_time = 0, total_write_time = 0, total_mat_time = 0, total_ref_time = 0,
         total_deref_time = 0, cg_partial;
    uint32_t total_cg_it = 0;
    struct stop_watch solver_time, ode_time, cg_time, part_solver, part_mat, write_time, ref_time, deref_time;

    init_stop_watch (&solver_time);
    init_stop_watch (&ode_time);
    init_stop_watch (&cg_time);
    init_stop_watch (&part_solver);
    init_stop_watch (&part_mat);
    init_stop_watch (&write_time);
    init_stop_watch (&ref_time);
    init_stop_watch (&deref_time);

    print_to_stdout_and_file ("Assembling Monodomain Matrix Begin\n");
    start_stop_watch (&part_mat);
    set_initial_conditions (the_monodomain_solver, the_grid, initial_v);
    set_discretization_matrix (the_monodomain_solver, the_grid);
    total_mat_time = stop_stop_watch (&part_mat);
    print_to_stdout_and_file ("Assembling Monodomain Matrix End\n");

    start_stop_watch (&solver_time);

    double cur_time = 0.0;

    int print_rate = configs->print_rate;

    bool abort_on_no_activity = the_monodomain_solver->abort_on_no_activity;
    double cg_error;
    uint32_t cg_iterations;

    set_spatial_stim (the_grid, stimuli_configs);

    if(has_extra_data)
        set_ode_extra_data(extra_data_config, the_grid, the_ode_solver);

    uint32_vector *refined_this_step = the_grid->refined_this_step;

    while (cur_time < finalT) {

        if (save_to_file) {

            if (count % print_rate == 0) {
                start_stop_watch (&write_time);

                sds tmp = sdsnew(configs->out_dir_name);
                sds c = sdsfromlonglong (count);
                tmp = sdscat (tmp, "/V_t_");
                tmp = sdscat (tmp, c);
                FILE *f1 = fopen (tmp, "w");
                activity = print_grid_and_check_for_activity (the_grid, f1, count);
                fclose (f1);
                sdsfree (tmp);
                sdsfree (c);

                total_write_time += stop_stop_watch (&write_time);

                if (abort_on_no_activity) {
                    if (!activity) {
                        print_to_stdout_and_file ("No activity, aborting simulation\n");
                        break;
                    }
                }
            }
        }

        count++;

        if (cur_time > 0.0) {
            update_ode_state_vector (the_ode_solver, the_grid, original_num_cells);
        }

        start_stop_watch (&ode_time);

        solve_all_volumes_odes(the_ode_solver, the_grid->num_active_cells, cur_time, ode_step, stimuli_configs);

        update_monodomain (original_num_cells, the_grid->num_active_cells, the_grid->active_cells, beta, cm, dt_edp,
                           the_ode_solver->sv, the_ode_solver->model_data.number_of_ode_equations, gpu);

        ode_total_time += stop_stop_watch (&ode_time);

        start_stop_watch (&cg_time);

        cg_iterations = conjugate_gradient (the_grid, max_its, cg_tol, jacobi, &cg_error);

        cg_partial = stop_stop_watch (&cg_time);

        cg_total_time += cg_partial;

        total_cg_it += cg_iterations;

        if (count % print_rate == 0) {
            print_to_stdout_and_file ("t = %lf, Iterations = "
                    "%" PRIu32 ", Error Norm = %e, Number of Cells:"
                    "%" PRIu32 ", Iterations time: %ld μs\n",
                    cur_time, cg_iterations, cg_error, the_grid->num_active_cells, cg_partial);
        }

        cur_time += dt_edp;

        if (adaptive) {
            redo_matrix = false;
            if (cur_time >= start_adpt_at) {

                if (count % refine_each == 0) {
                    start_stop_watch (&ref_time);
                    redo_matrix = refine_grid_with_bound (the_grid, refinement_bound, start_h);
                    total_ref_time += stop_stop_watch (&ref_time);
                }

                if (count % derefine_each == 0) {
                    start_stop_watch (&deref_time);
                    redo_matrix |= derefine_grid_with_bound (the_grid, derefinement_bound, max_h);
                    total_deref_time += stop_stop_watch (&deref_time);
                }
            }

            if (redo_matrix) {
                order_grid_cells (the_grid);

                set_spatial_stim (the_grid, stimuli_configs);
                if(has_extra_data)
                    set_ode_extra_data(extra_data_config, the_grid, the_ode_solver);

                update_cells_to_solve (the_grid, the_ode_solver);

                if (uint32_vector_size (refined_this_step) > 0) {
                    update_state_vectors_after_refinement (the_ode_solver, refined_this_step);
                }

                start_stop_watch (&part_mat);

                set_discretization_matrix (the_monodomain_solver, the_grid);

                total_mat_time += stop_stop_watch (&part_mat);
            }
        }
    }

    print_to_stdout_and_file ("Resolution Time: %ld μs\n", stop_stop_watch (&solver_time));
    print_to_stdout_and_file ("ODE Total Time: %ld μs\n", ode_total_time);
    print_to_stdout_and_file ("CG Total Time: %ld μs\n", cg_total_time);
    print_to_stdout_and_file ("Mat time: %ld μs\n", total_mat_time);
    print_to_stdout_and_file ("Refine time: %ld μs\n", total_ref_time);
    print_to_stdout_and_file ("Derefine time: %ld μs\n", total_deref_time);
    print_to_stdout_and_file ("Write time: %ld μs\n", total_write_time);
    print_to_stdout_and_file ("CG Total Iterations: %u\n", total_cg_it);
}


void set_spatial_stim (struct grid *the_grid, struct stim_config_hash *stim_configs) {

    struct stim_config *tmp = NULL;
    int n_active = the_grid->num_active_cells;

    for (int i = 0; i < stim_configs->size; i++) {
        for (struct stim_config_elt *e = stim_configs->table[i % stim_configs->size]; e != 0; e = e->next) {
            tmp = e->value;
            free(tmp->spatial_stim_currents);
            tmp->spatial_stim_currents = (Real*)malloc(sizeof(Real)*n_active);
            tmp->set_spatial_stim_fn(the_grid, tmp->stim_current, tmp->spatial_stim_currents, tmp->config_data.config);
        }
    }

}

void set_ode_extra_data(struct extra_data_config *config, struct grid *the_grid, struct ode_solver *the_ode_solver) {

    free(the_ode_solver->edo_extra_data);
    the_ode_solver->edo_extra_data = config->set_extra_data_fn(the_grid, config->config_data.config, &(the_ode_solver->extra_data_size));

}

void update_ode_state_vector (struct ode_solver *the_ode_solver, struct grid *the_grid, uint32_t max_number_of_cells) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    int n_edos = the_ode_solver->model_data.number_of_ode_equations;

    Real *sv = the_ode_solver->sv;

    if (the_ode_solver->gpu) {
#ifdef COMPILE_CUDA
        Real *vms;
        size_t mem_size = max_number_of_cells * sizeof (Real);

        vms = (Real *)malloc (mem_size);
        check_cuda_errors(cudaMemcpy (vms, sv, mem_size, cudaMemcpyDeviceToHost));

#pragma omp parallel for
        for (int i = 0; i < n_active; i++) {
            vms[ac[i]->sv_position] = (Real)ac[i]->v;
        }

        check_cuda_errors(cudaMemcpy (sv, vms, mem_size, cudaMemcpyHostToDevice));
        free (vms);
#endif
    } else {
#pragma omp parallel for
        for (uint32_t i = 0; i < n_active; i++) {
            sv[ac[i]->sv_position * n_edos] = (Real)ac[i]->v;
        }
    }
}

void save_old_cell_positions (struct grid *the_grid) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

#pragma omp parallel for
    for (uint32_t i = 0; i < n_active; i++) {
        ac[i]->sv_position = ac[i]->grid_position;
    }
}

void update_cells_to_solve (struct grid *the_grid, struct ode_solver *solver) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    if (solver->cells_to_solve) {
        free (solver->cells_to_solve);
    }

    solver->cells_to_solve = (uint32_t *)malloc (the_grid->num_active_cells * sizeof (uint32_t));
    uint32_t *cts = solver->cells_to_solve;

#pragma omp parallel for
    for (uint32_t i = 0; i < n_active; i++) {
        cts[i] = ac[i]->sv_position;
    }
}

void set_initial_conditions (struct monodomain_solver *the_solver, struct grid *the_grid, double initial_v) {

    double alpha, h;
    struct cell_node **ac = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;
    double beta = the_solver->beta;
    double cm = the_solver->cm;
    double dt = the_solver->dt;

#pragma omp parallel for private(alpha, h)
    for (int i = 0; i < active_cells; i++) {
        h = ac[i]->face_length;
        alpha = ALPHA (beta, cm, dt, h);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v * alpha;
    }
}

void set_discretization_matrix (struct monodomain_solver *the_solver, struct grid *the_grid) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    uint8_t max_elements = the_grid->num_cell_neighbours;

    initialize_diagonal_elements (the_solver, the_grid);

#pragma omp parallel for
    for (int i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements (the_solver, ac[i], ac[i]->south, 's', max_elements);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements (the_solver, ac[i], ac[i]->north, 'n', max_elements);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements (the_solver, ac[i], ac[i]->east, 'e', max_elements);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements (the_solver, ac[i], ac[i]->west, 'w', max_elements);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements (the_solver, ac[i], ac[i]->front, 'f', max_elements);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements (the_solver, ac[i], ac[i]->back, 'b', max_elements);
    }
}

void initialize_diagonal_elements (struct monodomain_solver *the_solver, struct grid *the_grid) {

    double alpha, h;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    double beta = the_solver->beta;
    double cm = the_solver->cm;
    double dt = the_solver->dt;

    uint8_t max_elements = the_grid->num_cell_neighbours;

#pragma omp parallel for private(alpha, h)
    for (int i = 0; i < num_active_cells; i++) {
        h = ac[i]->face_length;
        alpha = ALPHA (beta, cm, dt, h);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if (ac[i]->elements != NULL) {
            free (ac[i]->elements);
        }

        ac[i]->elements = new_element_array (max_elements);

        ac[i]->elements[0] = element;
    }
}

void fill_discretization_matrix_elements (struct monodomain_solver *the_solver, struct cell_node *grid_cell,
                                          void *neighbour_grid_cell, char direction, uint8_t max_elements) {

    uint32_t position;
    bool has_found;
    double h;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    double sigmaX = the_solver->sigma_x;
    double sigmaY = the_solver->sigma_y;
    double sigmaZ = the_solver->sigma_z;

    double sigmaX1 = (2.0f * sigmaX * sigmaX) / (sigmaX + sigmaX);
    double sigmaX2 = (2.0f * sigmaX * sigmaX) / (sigmaX + sigmaX);
    double sigmaY1 = (2.0f * sigmaY * sigmaY) / (sigmaY + sigmaY);
    double sigmaY2 = (2.0f * sigmaY * sigmaY) / (sigmaY + sigmaY);
    double sigmaZ1 = (2.0f * sigmaZ * sigmaZ) / (sigmaZ + sigmaZ);
    double sigmaZ2 = (2.0f * sigmaZ * sigmaZ) / (sigmaZ + sigmaZ);

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if (neighbour_grid_cell_level > grid_cell->cell_data.level) {
        if ((neighbour_grid_cell_type == 'w')) {
            has_found = false;
            while (!has_found) {
                if (neighbour_grid_cell_type == 'w') {
                    white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                    if (white_neighbor_cell->single_connector == NULL) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    }
    else {
        if (neighbour_grid_cell_level <= grid_cell->cell_data.level && (neighbour_grid_cell_type == 'w')) {
            has_found = false;
            while (!has_found) {
                if (neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
                    white_neighbor_cell = (struct transition_node *)(neighbour_grid_cell);
                    if (white_neighbor_cell->single_connector == 0) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    }

    // Tratamos somente os pontos interiores da malha.
    if (neighbour_grid_cell_type == CELL_NODE_TYPE) {

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if (black_neighbor_cell->active) {

            if (black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) {
                h = black_neighbor_cell->face_length;
            } else {
                h = grid_cell->face_length;
            }

            struct element *cell_elements = grid_cell->elements;
            position = black_neighbor_cell->grid_position;

            lock_cell_node (grid_cell);

            int el_counter = 1;

            while (el_counter < max_elements && cell_elements[el_counter].cell != NULL &&
                   cell_elements[el_counter].column != position) {
                el_counter++;
            }

            // TODO: maybe each element can have a different sigma!
            if (el_counter < max_elements && cell_elements[el_counter].cell == NULL) {

                struct element new_element;
                new_element.column = position;
                if (direction == 'n') { // Z direction
                    new_element.value = -sigmaZ1 * h;
                    cell_elements[0].value += (sigmaZ1 * h);
                } else if (direction == 's') { // Z direction
                    new_element.value = -sigmaZ2 * h;
                    cell_elements[0].value += sigmaZ2 * h;
                } else if (direction == 'e') { // Y direction
                    new_element.value = -sigmaY1 * h;
                    cell_elements[0].value += (sigmaY1 * h);
                } else if (direction == 'w') { // Y direction
                    new_element.value = -sigmaY2 * h;
                    cell_elements[0].value += (sigmaY2 * h);
                } else if (direction == 'f') { // X direction
                    new_element.value = -sigmaX1 * h;
                    cell_elements[0].value += (sigmaX1 * h);
                } else if (direction == 'b') { // X direction
                    new_element.value = -sigmaX2 * h;
                    cell_elements[0].value += (sigmaX2 * h);
                }

                new_element.cell = black_neighbor_cell;
                cell_elements[el_counter] = new_element;
            }
            unlock_cell_node (grid_cell);


            cell_elements = black_neighbor_cell->elements;
            position = grid_cell->grid_position;

            lock_cell_node (black_neighbor_cell);

            el_counter = 1;

            while (el_counter < max_elements && cell_elements[el_counter].cell != NULL &&
                   cell_elements[el_counter].column != position) {
                el_counter++;
            }

            if (el_counter < max_elements && cell_elements[el_counter].cell == NULL) {

                struct element new_element;
                new_element.column = position;
                if (direction == 'n') { // Z direction
                    new_element.value = -sigmaZ1 * h;
                    cell_elements[0].value += (sigmaZ1 * h);
                } else if (direction == 's') { // Z direction
                    new_element.value = -sigmaZ2 * h;
                    cell_elements[0].value += (sigmaZ2 * h);
                } else if (direction == 'e') { // Y direction
                    new_element.value = -sigmaY1 * h;
                    cell_elements[0].value += (sigmaY1 * h);
                } else if (direction == 'w') { // Y direction
                    new_element.value = -sigmaY2 * h;
                    cell_elements[0].value += (sigmaY2 * h);
                } else if (direction == 'f') { // X direction
                    new_element.value = -sigmaX1 * h;
                    cell_elements[0].value += (sigmaX1 * h);
                } else if (direction == 'b') { // X direction
                    new_element.value = -sigmaX2 * h;
                    cell_elements[0].value += (sigmaX2 * h);
                }

                new_element.cell = grid_cell;
                cell_elements[el_counter] = new_element;
            }

            unlock_cell_node (black_neighbor_cell);
        }
    }
}

void update_monodomain (uint32_t initial_number_of_cells, uint32_t num_active_cells, struct cell_node **active_cells,
                        double beta, double cm, double dt_edp, Real *sv, int n_equations_cell_model, bool use_gpu) {

    double h, alpha;

#ifdef COMPILE_CUDA
    Real *vms = NULL;
    size_t mem_size = initial_number_of_cells * sizeof (Real);

    if (use_gpu) {
        vms = (Real *)malloc (mem_size);
        check_cuda_errors(cudaMemcpy (vms, sv, mem_size, cudaMemcpyDeviceToHost));
    }
#endif

#pragma omp parallel for private(h, alpha)
    for (int i = 0; i < num_active_cells; i++) {
        h = active_cells[i]->face_length;
        alpha = ALPHA (beta, cm, dt_edp, h);

        if (use_gpu) {
#ifdef COMPILE_CUDA
            active_cells[i]->b = vms[active_cells[i]->sv_position] * alpha;
#endif
        } else {
            active_cells[i]->b = sv[active_cells[i]->sv_position * n_equations_cell_model] * alpha;
        }
    }
#ifdef COMPILE_CUDA
    free (vms);
#endif
}

//TODO: we have to pass more info to this function to print more information and print to a file as well
void print_solver_info(struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver,
                       struct grid *the_grid, struct user_options *options) {
    print_to_stdout_and_file ("System parameters: \n");
#if defined(_OPENMP)
    print_to_stdout_and_file ("Using OpenMP with %d threads\n", omp_get_max_threads ());
#endif
    if (the_ode_solver->gpu) {
        print_to_stdout_and_file ("Using GPU to solve ODEs\n");
    }

    print_to_stdout_and_file ("Initial V: %lf\n", the_ode_solver->model_data.initial_v);
    print_to_stdout_and_file ("Number of ODEs in cell model: %d\n", the_ode_solver->model_data.number_of_ode_equations);

    print_to_stdout_and_file ("Sigma X = %.10lf, Sigma Y = %.10lf, Sigma Z = %.10lf\n",
            the_monodomain_solver->sigma_x,
            the_monodomain_solver->sigma_y,
            the_monodomain_solver->sigma_z);

    print_to_stdout_and_file ("Initial N. of Elements = "
                    "%" PRIu32 "\n",
            the_grid->num_active_cells);
    print_to_stdout_and_file ("PDE time step = %lf\n", the_monodomain_solver->dt);
    print_to_stdout_and_file ("ODE min time step = %lf\n", the_ode_solver->min_dt);
    print_to_stdout_and_file ("Simulation Final Time = %lf\n", the_monodomain_solver->final_time);
    print_to_stdout_and_file ("Maximum CG iterations = %d\n", the_monodomain_solver->max_iterations);
    print_to_stdout_and_file ("CG tolerance = %e\n", the_monodomain_solver->tolerance);
    if (the_monodomain_solver->use_jacobi) {
        print_to_stdout_and_file ("Using Jacobi preconditioner\n");
    }
    if (the_grid->adaptive) {
        print_to_stdout_and_file ("Using adaptativity\n");
        print_to_stdout_and_file ("Refinement Bound = %lf\n", the_monodomain_solver->refinement_bound);
        print_to_stdout_and_file ("Derefinement Bound = %lf\n", the_monodomain_solver->derefinement_bound);
        print_to_stdout_and_file ("Refining each %d time steps\n", the_monodomain_solver->refine_each);
        print_to_stdout_and_file ("Derefining each %d time steps\n", the_monodomain_solver->derefine_each);
    }

    print_to_stdout_and_file ("Print Rate = %d\n", options->print_rate);

    if (options->out_dir_name != NULL) {
        print_to_stdout_and_file ("Saving to plain text output in %s dir\n", options->out_dir_name);
    } else {
        print_to_stdout_and_file ("The solution will not be saved\n");
    }


    if(options->stim_configs) {
        print_to_stdout_and_file("======================================================================\n");

        if (options->stim_configs->size == 1)
            print_to_stdout_and_file("Stimulus configuration:\n");
        else {
            print_to_stdout_and_file("Stimuli configuration:\n");
        }


        for (int i = 0; i < options->stim_configs->size; i++) {
            for (struct stim_config_elt *e = options->stim_configs->table[i % options->stim_configs->size];
                 e != 0; e = e->next) {

                print_to_stdout_and_file("Stimulus name: %s\n", e->key);
                print_to_stdout_and_file("Stimulus start: %lf\n", e->value->stim_start);
                print_to_stdout_and_file("Stimulus duration: %lf\n", e->value->stim_duration);
                print_to_stdout_and_file("Stimulus current: %lf\n", e->value->stim_current);
                print_to_stdout_and_file("Stimulus library: %s\n", e->value->config_data.library_file_path);
                print_to_stdout_and_file("Stimulus function: %s\n", e->value->config_data.function_name);
                struct string_hash *tmp = e->value->config_data.config;
                if (tmp->n == 1) {
                    print_to_stdout_and_file("Stimulus extra parameter:\n");
                } else if (tmp->n > 1) {
                    print_to_stdout_and_file("Stimulus extra parameters:\n");

                }

                STRING_HASH_PRINT_KEY_VALUE_LOG(tmp);

                print_to_stdout_and_file("======================================================================\n");

            }
        }
    }

    print_to_stdout_and_file("Domain configuration:\n");
    print_to_stdout_and_file("Domain name: %s\n", options->domain_config->domain_name);
    print_to_stdout_and_file("Domain initial Space Discretization: %lf um\n", options->domain_config->start_h);

    if (the_grid->adaptive) {
        print_to_stdout_and_file("Domain maximum Space Discretization: %lf um\n", options->domain_config->max_h);
        print_to_stdout_and_file("The adaptivity will start in time: %lf ms\n", the_monodomain_solver->start_adapting_at);
    }

    if(options->domain_config->config_data.config->n == 1) {
        print_to_stdout_and_file("Domain extra parameter:\n");
    }
    else if(options->domain_config->config_data.config->n > 1) {
        print_to_stdout_and_file("Domain extra parameters:\n");

    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(options->domain_config->config_data.config);
    print_to_stdout_and_file("======================================================================\n");

    if(options->extra_data_config) {
        print_to_stdout_and_file("Extra data ODE function configuration:\n");

        print_to_stdout_and_file("Extra data library: %s\n",options->extra_data_config->config_data.library_file_path);
        print_to_stdout_and_file("Extra data function: %s\n", options->extra_data_config->config_data.library_file_path);

        if (options->domain_config->config_data.config->n == 1) {
            print_to_stdout_and_file("Extra data parameter:\n");
        } else if (options->domain_config->config_data.config->n > 1) {
            print_to_stdout_and_file("Extra data parameters:\n");

        }

        STRING_HASH_PRINT_KEY_VALUE_LOG(options->extra_data_config->config_data.config);
        print_to_stdout_and_file("======================================================================\n");
    }


}

void configure_monodomain_solver_from_options(struct monodomain_solver *the_monodomain_solver,
                                              struct user_options *options) {

    assert(the_monodomain_solver);
    assert(options);

    the_monodomain_solver->tolerance = options->cg_tol;
    the_monodomain_solver->num_threads = options->num_threads;
    the_monodomain_solver->max_iterations = options->max_its;
    the_monodomain_solver->final_time = options->final_time;

    the_monodomain_solver->refine_each = options->refine_each;
    the_monodomain_solver->derefine_each = options->derefine_each;
    the_monodomain_solver->refinement_bound = options->ref_bound;
    the_monodomain_solver->derefinement_bound = options->deref_bound;

    the_monodomain_solver->abort_on_no_activity = options->abort_no_activity;

    the_monodomain_solver->dt = options->dt_edp;
    the_monodomain_solver->use_jacobi = options->use_jacobi;

    the_monodomain_solver->sigma_x = options->sigma_x;
    the_monodomain_solver->sigma_y = options->sigma_y;
    the_monodomain_solver->sigma_z = options->sigma_z;
    the_monodomain_solver->start_adapting_at = options->start_adapting_at;

}