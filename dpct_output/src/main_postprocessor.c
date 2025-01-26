#include <stdlib.h>

#include "3dparty/ini_parser/ini.h"
#include "config/config_parser.h"
#include "config/postprocessor_config.h"

static void print_usage(char **argv) {
	fprintf(stderr, "Usage: %s config_file.ini\n", argv[0]);
}

int main(int argc, char **argv) {

	if(argc != 2) {
		print_usage(argv);
		exit(EXIT_FAILURE);
	}

	char *config_file = argv[1];

	postprocess_list pp_functions = NULL;

	if(ini_parse(config_file, parse_preprocessor_config, (void*)&pp_functions) < 0) {
		fprintf(stderr, "Error: Can't load the config file %s\n", config_file);
		return EXIT_FAILURE;
	}

	for(int i = 0; i < arrlen(pp_functions); i++) {
		pp_functions[i]->function(pp_functions[i]->function_params);
	}

	return EXIT_SUCCESS;

}
