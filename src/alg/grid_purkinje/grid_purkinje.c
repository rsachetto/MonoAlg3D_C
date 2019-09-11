#include <assert.h>
#include <inttypes.h>
#include <float.h>

#include "grid_purkinje.h"

#include "../../single_file_libraries/stb_ds.h"
#include "../../utils/file_utils.h"

struct grid_purkinje *new_grid_purkinje () 
{
    struct grid_purkinje *result = (struct grid_purkinje *)malloc(sizeof(struct grid_purkinje));
    result->first_cell = NULL;
    result->purkinje_cells = NULL;

    result->num_active_purkinje_cells = result->number_of_purkinje_cells = 0;

    result->the_network = new_graph();

    return result;
}