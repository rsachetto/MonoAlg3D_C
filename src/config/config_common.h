//
// Created by sachetto on 14/10/17.
//

#ifndef MONOALG3D_CONFIG_COMMON_H
#define MONOALG3D_CONFIG_COMMON_H

#include "../common_types/common_types.h"
#include <stdbool.h>

struct config {
    void *handle;

    char *main_function_name;
    bool main_function_name_was_set;

    char *init_function_name;
    bool init_function_name_was_set;

    char *end_function_name;
    bool end_function_name_was_set;

    char *library_file_path;
    bool library_file_path_was_set;

    struct string_hash_entry *config_data;

    void *persistent_data;

    void *main_function;
    void *init_function;
    void *end_function;

    // used by save and restore state
    char **extra_function_names;
    void **extra_functions;

    real_cpu current_t;
};

#ifdef __cplusplus
extern "C" {
#endif
struct config *alloc_and_init_config_data();
void free_config_data(struct config *cm);
void init_config_functions(struct config *config, char *default_lib, char *config_type);
#ifdef __cplusplus
}
#endif

#define LOG_COMMON_CONFIG(tag, s)                                                                                                                              \
    do {                                                                                                                                                       \
        if((s) == NULL) {                                                                                                                                      \
            log_info(tag " no configuration.\n");                                                                                                              \
        } else {                                                                                                                                               \
            log_info(tag " configuration:\n");                                                                                                                 \
            log_info(tag " library = %s\n", (s)->library_file_path);                                                                                           \
            log_info(tag " main function = %s\n", (s)->main_function_name);                                                                                    \
                                                                                                                                                               \
            if((s)->init_function_name) {                                                                                                                      \
                log_info(tag " init function = %s\n", (s)->init_function_name);                                                                                \
            }                                                                                                                                                  \
                                                                                                                                                               \
            if((s)->end_function_name) {                                                                                                                       \
                log_info(tag " end function = %s\n", (s)->end_function_name);                                                                                  \
            }                                                                                                                                                  \
                                                                                                                                                               \
            for(int __i = 0; __i < arrlen((s)->extra_function_names); __i++) {                                                                                 \
                log_info(tag " extra function %d = %s\n", __i + 1, (s)->extra_function_names[__i]);                                                            \
            }                                                                                                                                                  \
            STRING_HASH_PRINT_KEY_VALUE_LOG(tag, (s)->config_data);                                                                                            \
        }                                                                                                                                                      \
    } while(0)

#define CALL_EXTRA_FUNCTIONS(fn, ...)                                                                                                                          \
    do {                                                                                                                                                       \
        uint32_t num_extra_fn = arrlen(config->extra_functions);                                                                                               \
                                                                                                                                                               \
        for(int __i = 0; __i < num_extra_fn; __i++) {                                                                                                          \
            ((fn *)(config->extra_functions[__i]))(__VA_ARGS__);                                                                                               \
        }                                                                                                                                                      \
    } while(0)

#endif // MONOALG3D_CONFIG_COMMON_H
