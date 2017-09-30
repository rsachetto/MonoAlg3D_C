#include <stdio.h>
#include <unitypes.h>
#include "grid.h"

int main() {

    struct grid *the_grid;

    FILE *f1;


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

    free_grid(the_grid);

    return 0;
}