#ifndef MONOALG3D_EXTRA_DATA_LIB_H
#define MONOALG3D_EXTRA_DATA_LIB_H

#include "../hash/string_hash.h"
#include "../alg/grid/grid.h"


void * set_extra_data_for_fibrosis_sphere(struct grid *the_grid, struct string_hash *extra_data_config, size_t *extra_data_bytes);
void * set_extra_data_for_fibrosis_plain(struct grid *the_grid, struct string_hash *extra_data_config, size_t *extra_data_bytes);

#endif //MONOALG3D_GRID_DOMAIN_H