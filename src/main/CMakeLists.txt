set(SOURCE_FILES_SOLVER
        stim_config_hash.c stim_config_hash.h
        stim_config.h stim_config.c
        ode_solver.c ode_solver.h
        monodomain_solver.h monodomain_solver.c
        linear_system_solver.c linear_system_solver.h
        config_parser.c config_parser.h
        constants.h)

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})