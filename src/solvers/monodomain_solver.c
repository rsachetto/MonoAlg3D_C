//
// Created by sachetto on 03/10/17.
//

#include "monodomain_solver.h"
#include "../gpu/gpu_common.h"
#include <omp.h>
#include <time.h>
#include "constants.h"

void init_solver(struct monodomain_solver* the_solver) {
    the_solver->cells_to_solve = NULL;
    the_solver->beta = 0.14;
    the_solver->cm = 1.0;

    the_solver->sigma_y = 0.0001334/2.5;
    the_solver->sigma_x = the_solver->sigma_y;
    the_solver->sigma_z = the_solver->sigma_y;

}

void solve_monodomain( struct grid *the_grid,  struct monodomain_solver *the_monodomain_solver, struct edo_solver *the_edo_solver, struct output_info *output_info ) {

    int refine_each = the_monodomain_solver->refine_each;
    int derefine_each = the_monodomain_solver->derefine_each;

    double cg_tol = the_monodomain_solver->tolerance;

    int np = the_monodomain_solver->num_threads;

    if(np == 0) np = 1;
    omp_set_num_threads(np);

    double initialV;
    bool redoMatrix = false;
    float *sv = NULL;
    size_t pitch = 0;

    bool activity = false;

    int max_its;

    //TODO: I dont think we need this flag anymore
    bool parallel = (np != 0);
    bool gpu = the_edo_solver->gpu;
    int count = 0;

    int maximumRefinementLevel;
    int maximumDerefinementeLevel;

    double refinement_bound = the_grid->refinement_bound;
    double derefinement_bound = the_grid->derefinement_bound;

    bool adaptive = the_grid->adaptive;
    bool save_to_file = (output_info->file_base_name != NULL);

    double h;
    int initRef;
    order_grid_cells(the_grid);
    the_grid->original_num_cells = the_grid->num_active_cells;

    double dt_edp = the_monodomain_solver->dt;
    double dt_edo = the_edo_solver->min_dt;
    double finalT = the_monodomain_solver->final_time;


    //TODO: create a file for stimulus definition!!
    double stim_start = 0.0;
    double stim_dur = the_monodomain_solver->stim_dur;
    double i_stim;

    int edo_method = the_edo_solver->method;


    if (gpu) {
        init_cuda_device(the_edo_solver->gpu_id);

        if(adaptive) {
            save_old_cell_positions(the_grid);
            update_cells_to_solve(the_grid, the_monodomain_solver);
        }
    }

    printf("Setting ODE's initial conditions\n");
    set_edo_initial_conditions(the_edo_solver);


    /*
    //TODO: we should decide how the stimulus will be handled
     initialV = INITIAL_V;
    if(isinf(globalArgs.stim_cur)) {
        if(globalArgs.benchmark) {
            i_stim = -50.0;
        }
        else if(globalArgs.use_plain || globalArgs.use_plain_with_sphere || globalArgs.use_human) {
            i_stim = -38;
        }
        else {
            i_stim = initialV;
        }

    }
    else {
        i_stim = globalArgs.stim_cur;
    }
    */

    clock_t mat_init, mat_final;

    printf("Assembling Monodomain Matrix\n");
    mat_init = clock();

    set_initial_conditions(the_monodomain_solver, the_grid);
    parallelSetDiscretizationMatrixMonodomain( beta, cm, dt_edp );


    mat_final = clock() - mat_init;
    cout << "Creating Matrix time: " << (double) mat_final/(double)CLOCKS_PER_SEC << " s" << endl;


    if(globalArgs.max_its > 0) {
        max_its = globalArgs.max_its;
    }
    else {
        max_its = numberOfCells;
    }

    double maxStep = globalArgs.max_dt_edo; //max time step allowed

    double abstol = globalArgs.abs_tol, reltol=globalArgs.rel_tol;

    printf("System parameters: \n");
    if(parallel) {
        printf("Using OpenMP with %d threads\n", np);
    }
    if(gpu) {
        printf("Using GPU to solve ODEs\n");
    }

    printf("Time discretization: %lf\n", dt_edp);
    printf("Initial V: %lf\n", initialV);
    printf("Initial Space Discretization: %lf um\n", globalArgs.start_h);

    if(adaptive) {
        printf("Initial Refinement Level: %d\n", getRefinementLevel(sideLength, globalArgs.start_h));
        printf("Minimum Space Discretization: %lf um\n", globalArgs.min_h);
        printf("Minimum Refinement Level: %d\n", getRefinementLevel(sideLength, globalArgs.min_h));
        printf("Maximum Space Discretization: %lf um\n", globalArgs.max_h);
        printf("Maximum Refinement Level: %d\n", getRefinementLevel(sideLength, globalArgs.max_h));
    }

    printf("Width = %lf um, height = %lf um, Depth = %lf \n", sideLength, sideLength, sideLength);
    printf("Volumes X = %d, Volumes Y = %d, Volumes Z = %d\n", int(ceil(sideLength/h)), int(ceil(sideLength/h)), int(ceil(sideLength/h)) );
    printf("Initial N. of Elements = %lld\n", numberOfCells);
    printf("PDE time step = %lf\n", dt_edp);
    printf("ODE solver edo_method: %s\n", getODEMethod(edo_method).c_str());
    printf("ODE min time step = %lf\n", dt_edo);
    if ( (edo_method == EULER_METHOD_ADPT) ) {
        printf("ODE max time step = %lf\n", maxStep);
        printf("Absolute tolerance for edo = %lf\n", abstol);
        printf("Relative tolerance for edo = %lf\n", reltol);
    }

    printf("Simulation Final Time = %lf\n", finalT);
    printf("Print Rate = %d\n", print_rate);
    printf("Stimulus start = %lf\n", stim_start);
    printf("Stimulus duration = %lf\n", stim_dur);
    printf("Stimulus value = %lf\n", i_stim);
    printf("Maximum CG iterations = %d\n", max_its);
    printf("CG tolerance = %e\n", cg_tol);
    if(globalArgs.use_jacobi) {
        printf("Using Jacobi preconditioner\n");
    }
    if (adaptive) {
        printf("Using adaptativity\n");
        printf("Refinement Bound = %lf\n", refinement_bound);
        printf("Derefinement Bound = %lf\n", derefinement_bound);
        printf("Refining each %d time steps\n", refine_each);
        printf("Derefining each %d time steps\n", derefine_each);

    }
    if (save_to_file) {
        if(!dirExists(globalArgs.outDirName)) {
            printf("%s does not exist!\n", globalArgs.outDirName);
            exit(0);
        }

#ifdef VTK
        printf("Saving to VTK output in %s dir\n", globalArgs.outDirName);
#else
        printf("Saving to plain text output in %s dir\n", globalArgs.outDirName);
#endif
    }
    else {
        printf("The solution will not be saved\n");
    }

    long ode_total_time = 0, cg_total_time = 0;
    int total_cg_it = 0;
    StopWatch solver_time, ode_time, cg_time;

    long total_write_time = 0, total_mat_time = 0, total_ref_time = 0, total_deref_time = 0;
    StopWatch part_solver, part_mat, write_time, ref_time, deref_time;

    solver_time.start();

    int ode_step = 1;

    if(edo_method == EULER_METHOD) {
        if(dt_edp >= dt_edo) {
            ode_step =(int)(dt_edp/dt_edo);
            printf("Solving EDO %d times before solving PDE\n", ode_step);
        }
        else {
            printf("WARNING: EDO time step is greater than PDE time step. Adjusting to EDO time step: %lf\n", dt_edo);
            dt_edp = dt_edo;
        }
    }

    fflush(stdout);

    if (!gpu) {
        set_stim_odes(stim_start, stim_dur);
    }

    if (edo_method ==  EULER_METHOD_ADPT  && !gpu ) {
        int aux = 0;
        CellNode *grid_cell;
        grid_cell = firstCell;
        set_global_odes_params(globalArgs.max_dt_edo, globalArgs.abs_tol, globalArgs.rel_tol);
        while( grid_cell != 0 ) {
            set_ode_params(grid_cell->od, dt_edo);
            grid_cell  = grid_cell->next;
            aux++;
        }
    }

    for( double cur_time = 0.0; cur_time < finalT; cur_time += dt_edp ) {

        if(save_to_file) {

            if(count % print_rate == 0) {
                write_time.start();

                stringstream file_name2 (stringstream::in | stringstream::out);
                //file_name2 << globalArgs.outDirName << "/V_t_" << count;

#ifdef VTK
                file_name2 << globalArgs.outDirName << "/V_t_" << count << ".vtu";

                activity = saveToVTK(file_name2.str(), count);

#else
                file_name2 << globalArgs.outDirName << "/V_t_" << count;
                activity = printOrdered(file_name2.str(), count);
#endif
                total_write_time += write_time.stop();


                if(globalArgs.use_plain || globalArgs.use_human || globalArgs.use_plain_with_sphere) {
                    if(!activity) {
                        cout << "No activity, aborting simulation" << endl;
                        break;
                    }
                }

            }
        }

        count++;

        if (cur_time > 0.0) {
            if(gpu) {
                if(adaptive) {
                    updateSvGPUDxAdpt(originalNumCells, numActiveCells, activeCells.data(), sv);
                }
                else{
                    updateSvGPU(numActiveCells, activeCells.data(), sv);
                }
            }
            else {
                setNewSV();
            }
        }

        ode_time.start();

        //TODO: remover (espiral)
        //if(t > 2.0) {
        //    stim_start = second_stim;
        //}


        if(edo_method == EULER_METHOD) {
            if(parallel || gpu) {
                if(gpu) {
                    if(adaptive) {
                        solveODEsDxAdptGPU(sv, cur_time, beta, cm, dt_edo, dt_edp, i_stim, stim_start, stim_dur, pitch, ode_step);
                        updateMonodomainDxAdpt(originalNumCells, numActiveCells, activeCells.data(), beta, cm, dt_edp, sv);
                    }
                    else {
                        solveODEsGPU(sv, cur_time, beta, cm, dt_edo, dt_edp, i_stim, stim_start, stim_dur, pitch, ode_step);
                        updateMonodomain(numActiveCells, activeCells.data(), beta, cm, dt_edp, sv);
                    }
                }
                else {
                    parallelSolveODEs(cur_time, beta, cm, dt_edo, dt_edp, i_stim, ode_step);
                }
            }
            else {
                solveODEs(cur_time, beta, cm, dt_edo, dt_edp, i_stim, ode_step);
            }

        }
        else if (edo_method == EULER_METHOD_ADPT) {
            if(parallel || gpu) {
                if (gpu) {
                    if(adaptive) {
                        cout << "solveODEsDtDxAdptGPU start" << endl;
                        solveODEsDtDxAdptGPU(sv, cur_time, beta, cm, dt_edp, i_stim, stim_start, stim_dur, globalArgs.min_dt_edo, globalArgs.max_dt_edo, pitch);
                        cout << "solveODEsDtDxAdptGPU end" << endl;

                        updateMonodomainDxAdpt(originalNumCells, numActiveCells, activeCells.data(), beta, cm, dt_edp, sv);

                    }
                    else {
                        solveODEsDtAdptGPU(sv, cur_time, beta, cm, dt_edp, i_stim, stim_start, stim_dur, globalArgs.min_dt_edo, globalArgs.max_dt_edo, pitch);
                        updateMonodomain(numActiveCells, activeCells.data(), beta, cm, dt_edp, sv);
                    }

                }
                else {
                    parallelSolveODEsAdpt(cur_time, beta, cm, dt_edp, i_stim);
                }
            }
            else {
                solveODEsAdpt(cur_time, beta, cm, dt_edp, i_stim);
            }

        }
        ode_total_time += ode_time.stop();
        cg_time.start();

        if(parallel) {
            parallelConjugateGradient(max_its, cg_tol, globalArgs.use_jacobi);
        }
        else {
            conjugateGradient(max_its, cg_tol, globalArgs.use_jacobi);
        }
        cg_total_time += cg_time.stop();

        total_cg_it += numberOfIterations;
        if(count % print_rate == 0) {
            printf("t = %lf, Iterations = %d, Error Norm = %e, Number of Cells: %d\n", cur_time, numberOfIterations, pError, numActiveCells);
        }


        if (adaptive) {
            redoMatrix = false;
            if ( cur_time > (stim_dur/2.0) ) {

                if(count % refine_each == 0) {
                    ref_time.start();
                    redoMatrix = refine(  globalArgs.min_h, refinement_bound );
                    total_ref_time += ref_time.stop();
                }


                if(count % derefine_each == 0) {
                    deref_time.start();
                    redoMatrix |= derefine( derefinement_bound,  globalArgs.max_h);
                    total_deref_time += deref_time.stop();
                }
            }

            if(redoMatrix) {
                ordergrid_cells(parallel || gpu);

                if(gpu) {
                    updateCellsToSolve();
                    if(refinedThisStep.size() > 0) {
                        updateAfterRefinement(sv, refinedThisStep.size(), &refinedThisStep[0], pitch, edo_method);
                    }
                }
                part_mat.start();
                if(parallel) {
                    parallelSetDiscretizationMatrixMonodomain( beta, cm, dt_edp );
                }
                else {
                    setDiscretizationMatrixMonodomain( beta, cm, dt_edp );
                }
                total_mat_time += part_mat.stop();
            }
        }

    }

    printf("Resolution Time: %ld us\n", solver_time.stop());;
    printf("ODE Total Time: %ld us\n", ode_total_time);
    printf("CG Total Time: %ld us\n", cg_total_time);
    cout << "Mat time: " << total_mat_time << " us" << endl;
    cout << "Refine time: " << total_ref_time << " us" << endl;
    cout << "Derefine time: " << total_deref_time << " us" << endl;
    cout << "Write time: " << total_write_time << " us" << endl;

    printf("CG Total Iterations: %d\n", total_cg_it);

}

void save_old_cell_positions(struct grid* the_grid) {

    #pragma omp parallel for
    for (int i = 0; i < the_grid->num_active_cells; i++) {
        the_grid->active_cells[i]->gpu_sv_position = the_grid->active_cells[i]->grid_position;
    }
}

void update_cells_to_solve(struct grid* the_grid, struct monodomain_solver *solver) {

    uint64_t n_active = the_grid->num_active_cells;

    if(solver->cells_to_solve) {
        free(solver->cells_to_solve);
    }

    solver->cells_to_solve = (uint64_t *)malloc(the_grid->num_active_cells*sizeof(uint64_t));

    #pragma omp parallel for
    for (uint64_t i = 0; i < n_active; i++) {
        solver->cells_to_solve[i] = the_grid->active_cells[i]->gpu_sv_position;
    }
}

void set_initial_conditions(struct monodomain_solver *the_solver, struct grid *the_grid) {

    double alpha, h;
    struct cell_node **ac = the_grid->active_cells;
    uint64_t active_cells = the_grid->num_active_cells;
    double beta = the_solver->beta;
    double cm = the_solver->cm;
    double dt = the_solver->dt;
    double initial_v = the_solver->initial_V;


    #pragma omp parallel for private (alpha, h)
    for (int i = 0; i < active_cells; i++) {
        h = ac[i]->face_length;
        alpha = ALPHA(beta, cm, dt, h);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v*alpha;
    }
}

void set_discretization_matrix(struct monodomain_solver *the_solver, struct grid *the_grid){

    void *neighbour_grid_cell;
    char direction;
    uint64_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    #pragma omp parallel for private (neighbour_grid_cell, direction)
    for (int i = 0; i < num_active_cells; i++) {

        //Computes and designates the flux due to south cells.
        neighbour_grid_cell = ac[i]->south;
        direction = 's';
        fill_discretization_matrix_elements( the_solver, ac[i], neighbour_grid_cell, direction);

        //Computes and designates the flux due to north cells.
        neighbour_grid_cell =  ac[i]->north;
        direction = 'n';
        fill_discretization_matrix_elements( the_solver, ac[i], neighbour_grid_cell, direction);

        //Computes and designates the flux due to east cells.
        neighbour_grid_cell =  ac[i]->east;
        direction = 'e';
        fill_discretization_matrix_elements( the_solver, ac[i], neighbour_grid_cell, direction);

        //Computes and designates the flux due to west cells.
        neighbour_grid_cell =  ac[i]->west;
        direction = 'w';
        fill_discretization_matrix_elements( the_solver, ac[i], neighbour_grid_cell, direction);

        neighbour_grid_cell =  ac[i]->front;
        direction = 'f';
        fill_discretization_matrix_elements( the_solver, ac[i], neighbour_grid_cell, direction);

        neighbour_grid_cell =  ac[i]->back;
        direction = 'b';
        fill_discretization_matrix_elements( the_solver, ac[i], neighbour_grid_cell, direction);

    }
}

void initialize_diagonal_elements(struct monodomain_solver *the_solver, struct grid *the_grid) {

    double alpha, h;
    uint64_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    double beta = the_solver->beta;
    double cm = the_solver->cm;
    double dt = the_solver->dt;


    //TODO: we need to decide how we will handle this linked list
    //Element *element;

    #pragma omp parallel for private (element, alpha, h)
    for (int i = 0; i < num_active_cells; i++) {
        h = ac[i]->face_length;
        alpha = ALPHA(beta, cm, dt, h);

        //TODO: we need to decide how we will handle this linked list
        /*
        element = new Element;
        element->column = ac[i]->grid_position;
        element->cell = ac[i];
        element->value = alpha;
        if (ac[i]->firstElement) {
            Element *aux = ac[i]->firstElement;
            while(aux) {
                Element *temp = aux;
                aux = aux->next;
                delete temp;
            }
        }
        ac[i]->firstElement = element;
         */
    }
}

void fill_discretization_matrix_elements(struct monodomain_solver *the_solver, struct cell_node *grid_cell, void* neighbour_grid_cell, char direction) {

    int position;
    bool has_found;
    double h;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    double sigmaX = the_solver->sigma_x;
    double sigmaY = the_solver->sigma_y;
    double sigmaZ = the_solver->sigma_z;


    double sigmaX1 = (2.0 * sigmaX * sigmaX) / (sigmaX + sigmaX);
    double sigmaX2 = (2.0 * sigmaX * sigmaX) / (sigmaX + sigmaX);
    double sigmaY1 = (2.0 * sigmaY * sigmaY) / (sigmaY + sigmaY);
    double sigmaY2 = (2.0 * sigmaY * sigmaY) / (sigmaY + sigmaY);
    double sigmaZ1 = (2.0 * sigmaZ * sigmaZ) / (sigmaZ + sigmaZ);
    double sigmaZ2 = (2.0 * sigmaZ * sigmaZ) / (sigmaZ + sigmaZ);

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    //Acha uma célula real que está no caixo enviado como vizinho
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data*)(neighbour_grid_cell))->level;
    char neighbour_grid_cell_type = ((struct basic_cell_data*)(neighbour_grid_cell))->type;

    if (neighbour_grid_cell_level > grid_cell->cell_data.level ) {
        if((neighbour_grid_cell_type == 'w') ) {
            white_neighbor_cell = (struct transition_node*)neighbour_grid_cell;
            has_found = false;
            while( !has_found ) {
                if( neighbour_grid_cell_type == 'w' ) {
                    white_neighbor_cell = (struct transition_node*)neighbour_grid_cell;
                    if( white_neighbor_cell->single_connector == NULL ) {
                        has_found = true;
                    }
                    else {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data*)(neighbour_grid_cell))->type;
                    }
                }
                else {
                    break;
                }
            }

        }
    }
        //Aqui, a célula vizinha tem um nivel de refinamento menor, entao eh mais simples.
    else {
        if(neighbour_grid_cell_level <= grid_cell->cell_data.level && (neighbour_grid_cell_type == 'w') ) {
            white_neighbor_cell = (struct transition_node*)(neighbour_grid_cell);
            has_found = false;
            while( !has_found ) {
                if( neighbour_grid_cell_type == 'w' ) {
                    white_neighbor_cell = (struct transition_node*)(neighbour_grid_cell);
                    if( white_neighbor_cell->single_connector == 0 ) {
                        has_found = true;
                    }
                    else {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data*)(neighbour_grid_cell))->type;
                    }
                }
                else {
                    break;
                }
            }
        }
    }

    //Tratamos somente os pontos interiores da malha.
    if( neighbour_grid_cell_type == 'b'  ) {

        black_neighbor_cell = (struct cell_node *) (neighbour_grid_cell);
        
        if (black_neighbor_cell->active) {

            if (black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) {
                h = black_neighbor_cell->face_length;
            } else {
                h = grid_cell->face_length;
            }

            //cout << "Update cell " << grid_cell->gridPosition << " and " << black_neighbor_cell->gridPosition << "  " << direction << endl;

            //Descobrimos a coluna que temos que preencher com o vizinho
            Element *element;
            element = grid_cell->firstElement;
            position = black_neighbor_cell->gridPosition;

            grid_cell->lock();

            while (element != 0 && element->column != position) {
                element = element->next;
            }

            //TODO: Cada elemento pode ter um sigma diferente
            if (element == NULL) {

                element = new
                Element;
                element->column = position;
                if (direction == 'n') { //Z direction
                    element->value = -sigmaZ1 * h;
                    grid_cell->firstElement->value += (sigmaZ1 * h);
                } else if (direction == 's') { //Z direction
                    element->value = -sigmaZ2 * h;
                    grid_cell->firstElement->value += (sigmaZ2 * h);
                } else if (direction == 'e') { //Y direction
                    element->value = -sigmaY1 * h;
                    grid_cell->firstElement->value += (sigmaY1 * h);
                } else if (direction == 'w') { //Y direction
                    element->value = -sigmaY2 * h;
                    grid_cell->firstElement->value += (sigmaY2 * h);
                } else if (direction == 'f') { //X direction
                    element->value = -sigmaX1 * h;
                    grid_cell->firstElement->value += (sigmaX1 * h);
                } else if (direction == 'b') {//X direction
                    element->value = -sigmaX2 * h;
                    grid_cell->firstElement->value += (sigmaX2 * h);
                }

                element->cell = black_neighbor_cell;
                Element *auxiliarElement;
                auxiliarElement = grid_cell->firstElement->next;
                grid_cell->firstElement->next = element;
                element->next = auxiliarElement;

            }
            grid_cell->unlock();

            //preenchemos a outra parte (a matrix é simetrica)
            element = black_neighbor_cell->firstElement;
            position = grid_cell->gridPosition;


            black_neighbor_cell->lock();

            while (element != 0 && element->column != position) {
                element = element->next;
            }


            if (element == NULL) {

                element = new Element;
                element->column = position;

                if (direction == 'n') { //Z direction
                    element->value = -sigmaZ1 * h;
                    black_neighbor_cell->firstElement->value += (sigmaZ1 * h);
                } else if (direction == 's') { //Z direction
                    element->value = -sigmaZ2 * h;
                    black_neighbor_cell->firstElement->value += (sigmaZ2 * h);
                } else if (direction == 'e') { //Y direction
                    element->value = -sigmaY1 * h;
                    black_neighbor_cell->firstElement->value += (sigmaY1 * h);
                } else if (direction == 'w') { //Y direction
                    element->value = -sigmaY2 * h;
                    black_neighbor_cell->firstElement->value += (sigmaY2 * h);
                } else if (direction == 'f') { //X direction
                    element->value = -sigmaX1 * h;
                    black_neighbor_cell->firstElement->value += (sigmaX1 * h);
                } else if (direction == 'b') { //X direction
                    element->value = -sigmaX2 * h;
                    black_neighbor_cell->firstElement->value += (sigmaX2 * h);
                }

                //TODO: barreira
                element->cell = grid_cell;
                Element *auxiliarElement;
                auxiliarElement = black_neighbor_cell->firstElement->next;
                black_neighbor_cell->firstElement->next = element;
                element->next = auxiliarElement;

            }

            black_neighbor_cell->unlock();
        }
    }
}