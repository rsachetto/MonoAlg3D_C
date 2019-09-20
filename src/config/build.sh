CONFIG_SOURCE_FILES="config_common.c stim_config.c domain_config.c
              purkinje_config.c assembly_matrix_config.c extra_data_config.c
   	          linear_system_solver_config.c
		          config_parser.c	save_mesh_config.c
		          save_state_config.c	restore_state_config.c
		          update_monodomain_config.c"
		          
CONFIG_HEADER_FILES="config_common.h stim_config.h domain_config.h
              purkinje_config.h assembly_matrix_config.h extra_data_config.h
   	          linear_system_solver_config.h
		          config_parser.h	save_mesh_config.h
		          save_state_config.h	restore_state_config.h
		          update_monodomain_config.h"

COMPILE_STATIC_LIB "config" "$CONFIG_SOURCE_FILES" "$CONFIG_HEADER_FILES"
