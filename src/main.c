#include <stdio.h>
#include <unitypes.h>
#include "grid/grid/grid.h"
#include "solvers/ode_solver.h"
#include "solvers/monodomain_solver.h"

int main(int argc, char **argv) {

    struct grid *the_grid;
    struct monodomain_solver *edp_solver;
    struct ode_solver *ode_solver;
    struct output_utils *output_info;

    Real dt = 0.02f;
    double start_h = 500.0f;
    double max_h = 500.0f;
    double min_h = start_h;

    the_grid = (struct grid*)malloc(sizeof(struct grid));
    edp_solver = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));

    output_info = new_output_utils(100, NULL);
    ode_solver = new_ode_solver("/home/sachetto/MonoAlg3D_C/model_lib/libten_tusscher_2006_cpu.so");
    ode_solver->gpu = false;
    ode_solver->min_dt = dt;
    ode_solver->stim_duration = 2.0;
    ode_solver->stim_start = 0.0;
    ode_solver->method = EULER_METHOD;
    ode_solver->stim_current = -50.0f;

    init_solver(edp_solver);

    edp_solver->tolerance = 1e-16;
    edp_solver->use_jacobi = true;
    edp_solver->num_threads = 1;
    edp_solver->dt = dt;
    edp_solver->adaptive = false;
    edp_solver->start_h = start_h;
    edp_solver->min_h = min_h;
    edp_solver->max_h = max_h;
    edp_solver->final_time = 20.0;
    edp_solver->abort_on_no_activity = false;
    edp_solver->use_jacobi = true;

    edp_solver->max_iterations = 200;

    initialize_grid_with_benchmark_mesh(the_grid, start_h);

    solve_monodomain(the_grid, edp_solver, ode_solver, output_info);

    free_output_utils(output_info);
    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);

    free(edp_solver);

    return 0;
}

