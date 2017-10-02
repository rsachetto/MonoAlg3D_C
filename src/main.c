#include <stdio.h>
#include <unitypes.h>
#include "grid/grid.h"
#include "opts.h"

int main(int argc, char **argv) {

    struct grid *the_grid;

    FILE *f1;

    struct user_args user_args;

    parse_options(argc, argv, &user_args);


    the_grid = (struct grid*)malloc(sizeof(struct grid));

    initialize_grid_with_mouse_mesh(the_grid, "../../Dropbox/Universidade/Pesquisa/alg/ALG3D/MonoAlg3D_v2_circle_images/meshes/mouse.alg");

    f1 = fopen("V_t_2", "w");
    print_grid(the_grid, f1);
    fclose(f1);

    refine_grid(the_grid, 1);
    f1 = fopen("V_t_3", "w");
    print_grid(the_grid, f1);
    fclose(f1);


    free_grid(the_grid);

    return 0;
}