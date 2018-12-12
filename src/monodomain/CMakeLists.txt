set(SOURCE_FILES_SOLVER
        output_utils.h
        output_utils.c
        ode_solver.c ode_solver.h
        monodomain_solver.h monodomain_solver.c
        constants.h)

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})

IF(MSVC)
	TARGET_LINK_LIBRARIES(solvers getopt)
ELSE()
	TARGET_LINK_LIBRARIES(solvers alg)
ENDIF()