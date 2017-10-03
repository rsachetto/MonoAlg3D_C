#include <stdio.h>
#include <unitypes.h>
#include "grid/grid.h"
#include "solvers/edo_solver.h"
#include "solvers/monodomain_solver.h"
#include "utils/output_info.h"
#include "utils/opts.h"

int main(int argc, char **argv) {

    struct grid *the_grid;
    struct monodomain_solver *edp_solver;
    struct edo_solver *edo_solver;
    struct output_info *output_info;


    FILE *f1;

    struct user_args user_args;

    parse_options(argc, argv, &user_args);


    the_grid = (struct grid*)malloc(sizeof(struct grid));
    edp_solver = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));
    edo_solver = (struct edo_solver*)malloc(sizeof(struct edo_solver));
    output_info = (struct output_info*)malloc(sizeof(struct output_info));


    //initialize_grid_with_mouse_mesh(the_grid, "../../Dropbox/Universidade/Pesquisa/alg/ALG3D/MonoAlg3D_v2_circle_images/meshes/mouse.alg");
    //initialize_grid_with_rabbit_mesh(the_grid, "../../Dropbox/Universidade/Pesquisa/alg/ALG3D/MonoAlg3D_v2_circle_images/meshes/rabheart.alg");
    //initialize_grid_with_benchmark_mesh(the_grid, 100.0);
    //initialize_grid_with_plain_mesh(the_grid, 10000.0, 100.0, 1);
    //f1 = fopen("V_t_1", "w");
    //print_grid(the_grid, f1);
    //fclose(f1);

    //clean_grid(the_grid);

    //initialize_grid_with_plain_fibrotic_mesh(the_grid, 10000, 100, 1, 0.45);
    //f1 = fopen("V_t_2", "w");
    //print_grid(the_grid, f1);
    //fclose(f1);

    //clean_grid(the_grid);

    float plain_center = 20050.0f; // this center works using a 40000um mesh
    float sphere_radius = 14000.0f;

    float bz_radius = 16000.0f;
    float bz_size = bz_radius - sphere_radius;

    init_solver(edp_solver);

    //TODO: we should think in a file format to handle mesh information
    initialize_grid_with_plain_and_sphere_fibrotic_mesh (the_grid, 40000.0f, 100.0f, 2, 0.25, plain_center, sphere_radius, bz_size, bz_radius);

    //TODO: now we have to implement this function
    solve_monodomain(the_grid, edp_solver, edo_solver, output_info);

//    f1 = fopen("V_t_3", "w");
//    print_grid(the_grid, f1);
//    fclose(f1);
//
//    clean_grid(the_grid);
//    free(the_grid);

    return 0;
}

