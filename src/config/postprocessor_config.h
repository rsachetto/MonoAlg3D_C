//
// Created by sachetto on 13/10/17.
//

#ifndef POSTPROCESSOR_CONFIG_H
#define POSTPROCESSOR_CONFIG_H

#include "../common_types/common_types.h"
#include "config_common.h"

#define DEFAULT_POSTPROCESSOR_LIB "shared_libs/libdefault_postprocess.so"

#define POSTPROCESS(name) int name(struct config *config)
typedef POSTPROCESS(postprocess_fn);

struct postprocess_function {
    postprocess_fn *function;
    struct config *function_params;
};

typedef struct postprocess_function ** postprocess_list;

struct postprocess_function *new_postprocess_function();
void init_postprocess_function(struct postprocess_function *function, const char *f_name);

#endif //POSTPROCESSOR_CONFIG_H
