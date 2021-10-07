CONFIG_SOURCE_FILES="config_common.c config_parser.c postprocessor_config.c"
		          
CONFIG_HEADER_FILES="config_common.h stim_config.h domain_config.h
                     purkinje_config.h assembly_matrix_config.h extra_data_config.h
   	                 linear_system_solver_config.h
		             config_parser.h save_mesh_config.h
		             save_state_config.h restore_state_config.h
		             update_monodomain_config.h postprocessor_config.h"

COMPILE_STATIC_LIB "config" "$CONFIG_SOURCE_FILES" "$CONFIG_HEADER_FILES"
