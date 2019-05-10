set(SOURCE_FILES_CONFIG
		config_common.c config_common.h
		stim_config.h stim_config.c
		domain_config.c domain_config.h
		purkinje_config.c purkinje_config.h
		assembly_matrix_config.c assembly_matrix_config.h
		extra_data_config.c extra_data_config.h
		linear_system_solver_config.c linear_system_solver_config.h
		config_parser.c config_parser.h
		save_mesh_config.c save_mesh_config.h
		save_state_config.c save_state_config.h
		restore_state_config.c restore_state_config.h
		update_monodomain_config.c update_monodomain_config.h
        )

ADD_LIBRARY(config STATIC ${SOURCE_FILES_CONFIG})

TARGET_LINK_LIBRARIES(config string utils dl)
