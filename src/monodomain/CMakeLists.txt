set(SOURCE_FILES_MONODOMAIN_SOLVER
        monodomain_solver.h monodomain_solver.c
        constants.h)

ADD_LIBRARY(monodomain STATIC ${SOURCE_FILES_MONODOMAIN_SOLVER})

TARGET_LINK_LIBRARIES(monodomain alg)
