#TODO: check zlib??
VTK_UTILS_SOURCE_FILES="pvd_utils.c data_utils.c vtk_polydata_grid.c vtk_unstructured_grid.c"
VTK_UTILS_HEADER_FILES="pvd_utils.h data_utils.h vtk_polydata_grid.h vtk_unstructured_grid.h"

COMPILE_STATIC_LIB "vtk_utils" "$VTK_UTILS_SOURCE_FILES" "$VTK_UTILS_HEADER_FILES"