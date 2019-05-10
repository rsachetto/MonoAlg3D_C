set(SOURCE_FILES_ALG
        grid/grid.h
        grid/grid.c
        grid/grid_refinement.c
        grid/grid_derefinement.c
        cell/cell.h
        cell/cell.c
        cell/cell_derefinement.c
        cell/cell_refinement.c
        )

ADD_LIBRARY(alg STATIC ${SOURCE_FILES_ALG})
TARGET_LINK_LIBRARIES(alg graph)
set_property(TARGET alg PROPERTY POSITION_INDEPENDENT_CODE ON)
