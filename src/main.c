#include <stdio.h>
#include <unitypes.h>
#include "grid.h"
#include "opts.h"

int main(int argc, char **argv) {

    struct grid *the_grid;

    FILE *f1;

    struct user_args user_args;

    parse_options(argc, argv, &user_args);


    the_grid = (struct grid*)malloc(sizeof(struct grid));

    initialize_grid(the_grid, 12800);

    construct_grid(the_grid);
    f1 = fopen("V_t_1", "w");
    print_grid(the_grid, f1);
    fclose(f1);

    refine_grid(the_grid, 1);
    f1 = fopen("V_t_2", "w");
    print_grid(the_grid, f1);
    fclose(f1);

    refine_grid(the_grid, 1);
    f1 = fopen("V_t_3", "w");
    print_grid(the_grid, f1);
    fclose(f1);


    refine_grid_cell_at(the_grid, 0);
    f1 = fopen("V_t_4", "w");
    print_grid(the_grid, f1);
    fclose(f1);

    refine_grid_cell_at(the_grid, 0);
    f1 = fopen("V_t_5", "w");
    print_grid(the_grid, f1);
    fclose(f1);

    refine_grid_cell_at(the_grid, 1);
    f1 = fopen("V_t_6", "w");
    print_grid(the_grid, f1);
    fclose(f1);


    derefine_all_grid(the_grid);
    f1 = fopen("V_t_7", "w");
    print_grid(the_grid, f1);
    fclose(f1);

    derefine_all_grid(the_grid);
    f1 = fopen("V_t_8", "w");
    print_grid(the_grid, f1);
    fclose(f1);

    derefine_all_grid(the_grid);
    f1 = fopen("V_t_9", "w");
    print_grid(the_grid, f1);
    fclose(f1);


    free_grid(the_grid);

    return 0;
}