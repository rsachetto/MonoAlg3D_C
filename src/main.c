#include <stdio.h>
#include <unitypes.h>
#include "grid/grid.h"
#include "solvers/ode_solver.h"
#include "solvers/monodomain_solver.h"
#include "utils/output_utils.h"
#include "utils/opts.h"

int main(int argc, char **argv) {

    struct grid *the_grid;
    struct monodomain_solver *edp_solver;
    struct ode_solver *ode_solver;
    struct output_utils *output_info;

    FILE *f1;

    struct user_args user_args;

    parse_options(argc, argv, &user_args);

    the_grid = (struct grid*)malloc(sizeof(struct grid));
    edp_solver = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));

    output_info = (struct output_utils*)malloc(sizeof(struct output_utils));

    ode_solver = new_ode_solver("/home/sachetto/MonoAlg3D_C/model_lib/libbondarenko_2004_cpu.so");
    ode_solver->gpu = false;

    int n = 2;
    set_ode_initial_conditions_for_all_volumes(ode_solver, n);

    printf("Cell model size %d and initial V %lf\n", ode_solver->model_data.number_of_ode_equations, ode_solver->model_data.initial_v);


    for(int i = 0; i < ode_solver->model_data.number_of_ode_equations*n; i++) {
        printf("%lf\n", ode_solver->sv[i]);
    }


//    init_solver(edp_solver);
//    edp_solver->num_threads = 1;
//    edp_solver->dt = 0.02;
//
//    //TODO: this should be provided by a user provided function
//    initialize_grid_with_plain_mesh(the_grid, 1000.0, 100.0f, 3);
//
//    solve_monodomain(the_grid, edp_solver, ode_solver, output_info);


    return 0;
}

