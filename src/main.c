#include <stdio.h>
#include <unitypes.h>
#include "grid.h"

int main() {

    struct grid *the_grid;

    the_grid = (struct grid*)malloc(sizeof(struct grid));

    initialize_grid(the_grid);
    construct_grid(the_grid);

    print_grid(the_grid, stdout);
    free_grid(the_grid);

    return 0;
}