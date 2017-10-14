set(SOURCE_FILES_ALG
        grid/grid.h
        grid/grid.c
        grid/grid_refinement.c
        grid/grid_derefinement.c
        cell/cell.h
        cell/cell.c
        cell/cell_derefinement.c
        cell/cell_refinement.c)

ADD_LIBRARY(alg SHARED ${SOURCE_FILES_ALG})