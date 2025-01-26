#ifndef __COMMON_H
#define __COMMON_H 

#include "../config/config_parser.h"

struct user_options *load_options_from_file(char *config_file);
int run_simulation_with_config(struct user_options *options, char *out_dir);

#endif /* __COMMON_H */
