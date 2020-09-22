ALG_SOURCE_FILES="grid/grid.c grid/grid_refinement.c grid/grid_derefinement.c cell/cell.c cell/cell_derefinement.c cell/cell_refinement.c  grid_purkinje/grid_purkinje.c"
ALG_HEADER_FILES="grid/grid.h cell/cell.h grid_purkinje/grid_purkinje.h"

COMPILE_STATIC_LIB "alg" "$ALG_SOURCE_FILES" "$ALG_HEADER_FILES"
