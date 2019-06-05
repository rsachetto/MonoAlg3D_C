set(SOURCE_FILES_SOLVER
        ode_solver.c ode_solver.h
        monodomain_solver.h monodomain_solver.c
        constants.h)

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})

TARGET_LINK_LIBRARIES(solvers alg)
