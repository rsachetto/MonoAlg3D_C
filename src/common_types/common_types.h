#ifndef MONOALG3D_COMMON_TYPES_H
#define MONOALG3D_COMMON_TYPES_H

#include <stdint.h>
#include "../monodomain/constants.h"

struct line {
    uint64_t source;
    uint64_t destination;
};

struct point_3d {
    real_cpu x, y, z;
};
#define POINT3D(x,y,z) (struct point_3d) {x,y,z}
#define ZERO_POINT3D (struct point_3d) {0,0,0}

struct point_hash_entry {
    struct point_3d key;
    float value;
};

struct string_hash_entry {
    char *key;
    char *value;
};

struct point_voidp_hash_entry {
    struct point_3d key;
    void *value;
};

struct string_voidp_hash_entry {
    char *key;
    void *value;
};

typedef uint32_t * ui32_array;
typedef  struct element * element_array;
typedef float * f32_array;
typedef struct point_3d * point3d_array;
typedef struct line * line_array;
typedef int64_t * int64_array;
typedef int * int_array;
typedef char** string_array;

struct vtk_files {
    string_array files_list;
    f32_array    timesteps;
};

#define STRING_HASH_PRINT_KEY_VALUE(d)                                                                                 \
    do {                                                                                                               \
        for(long i = 0; i < shlen(d); i++) {                                                                            \
            struct string_hash_entry e = d[i];                                                                         \
            printf("%s = %s\n", e.key, e.value);                                                                       \
        }                                                                                                              \
    } while(0)

#define STRING_HASH_PRINT_KEY_VALUE_LOG(d)                                                                             \
    do {                                                                                                               \
        for(long i = 0; i < shlen(d); i++) {                                                                            \
            struct string_hash_entry e = d[i];                                                                         \
            print_to_stdout_and_file("%s = %s\n", e.key, e.value);                                                     \
        }                                                                                                              \
    } while(0)

#define STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE(d, fn)                                                         \
    do {                                                                                                               \
        for(long i = 0; i < hmlen(d); i++) {                                                                            \
            struct string_voidp_hash_entry e = d[i];                                                                   \
            fn(e.value);                                                                                               \
        }                                                                                                              \
    } while(0)                                                                                                         \

#define STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE_AND_KEY(d, fn)                                                 \
    do {                                                                                                               \
        for(long i = 0; i < hmlen(d); i++) {                                                                            \
            struct string_voidp_hash_entry e = d[i];                                                                   \
            fn(e.value, e.key);                                                                                        \
        }                                                                                                              \
    }                                                                                                                  \
    while(0)

#define STIM_CONFIG_HASH_FOR_INIT_FUNCTIONS(d)                                                                         \
    do {                                                                                                               \
        for(long i = 0; i < hmlen(d); i++) {                                                                           \
            struct string_voidp_hash_entry e = d[i];                                                                   \
            init_config_functions(e.value, "./shared_libs/libdefault_stimuli.so", e.key);                              \
        }                                                                                                              \
    }                                                                                                                  \
    while(0)

#endif