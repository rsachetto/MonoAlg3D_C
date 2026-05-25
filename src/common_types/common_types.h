#ifndef MONOALG3D_COMMON_TYPES_H
#define MONOALG3D_COMMON_TYPES_H

#include "../3dparty/sds/sds.h"
#include <stdbool.h>
#include <stdint.h>

#define Pragma(x) _Pragma(#x)
#define OMP(directive) Pragma(omp directive)

typedef double real_cpu;

// Precision to be used for the calculations on the GPU
#ifdef CELL_MODEL_REAL_DOUBLE
typedef double real;
#else
typedef float real;
#endif

#define EPS (real)1e-16

#define ALLOCATE_3D_ARRAY(array, type, size_x, size_y, size_z)                                                                                                 \
    do {                                                                                                                                                       \
        array = (type ***)malloc(size_x * sizeof(type **));                                                                                                    \
                                                                                                                                                               \
        if(array == NULL) {                                                                                                                                    \
            log_error_and_exit("Memory allocation failed for first dimension\n");                                                                              \
        }                                                                                                                                                      \
                                                                                                                                                               \
        for(int i = 0; i < size_x; ++i) {                                                                                                                      \
            array[i] = (type **)malloc(size_y * sizeof(type *));                                                                                               \
            if(array[i] == NULL) {                                                                                                                             \
                log_error_and_exit("Memory allocation failed for second dimension\n");                                                                         \
            }                                                                                                                                                  \
                                                                                                                                                               \
            for(int j = 0; j < size_y; ++j) {                                                                                                                  \
                array[i][j] = (type *)calloc(size_z, sizeof(type));                                                                                            \
                if(array[i][j] == NULL) {                                                                                                                      \
                    log_error_and_exit("Memory allocation failed for third dimension\n");                                                                      \
                }                                                                                                                                              \
            }                                                                                                                                                  \
        }                                                                                                                                                      \
    } while(0)

#define FREE_3D_ARRAY(array, size_x, size_y)                                                                                                                   \
    do {                                                                                                                                                       \
        for(int i = 0; i < size_x; ++i) {                                                                                                                      \
            for(int j = 0; j < size_y; ++j) {                                                                                                                  \
                free(array[i][j]);                                                                                                                             \
            }                                                                                                                                                  \
            free(array[i]);                                                                                                                                    \
        }                                                                                                                                                      \
        free(array);                                                                                                                                           \
    } while(0)

#define MALLOC_BYTES(type, bytes) (type *)malloc(bytes)
#define MALLOC_ONE_TYPE(type) (type *)malloc(sizeof(type))
#define MALLOC_ARRAY_OF_TYPE(type, n) (type *)malloc(sizeof(type) * (n))

#define CALLOC_ONE_TYPE(type) (type *)calloc(1, sizeof(type))
#define CALLOC_ARRAY_OF_TYPE(type, n) (type *)calloc(n, sizeof(type))

#define ALLOCATE_MESH_INFO(grid_cell, mesh_info_struct)                                                                                                        \
    do {                                                                                                                                                       \
        size_t __size__ = sizeof(struct mesh_info_struct);                                                                                                     \
        (grid_cell)->mesh_extra_info = malloc(__size__);                                                                                                       \
        (grid_cell)->mesh_extra_info_size = __size__;                                                                                                          \
    } while(0)

#define MESH_INFO_DATA(grid_cell, mesh_info_struct, data_name) ((struct mesh_info_struct *)(grid_cell)->mesh_extra_info)->data_name

enum simulation_status {
    RESTART_SIMULATION,
    END_SIMULATION,
    SIMULATION_FINISHED,
};

struct time_info {
    real_cpu current_t;
    real_cpu final_t;
    real_cpu dt;
    int iteration;
};

#define TIME_INFO(it, ct, dt, lt)                                                                                                                              \
    (struct time_info) {                                                                                                                                       \
        it, ct, dt, lt                                                                                                                                         \
    }
#define ZERO_TIME_INFO TIME_INFO(0, 0, 0, 0)

struct line {
    uint64_t source;
    uint64_t destination;
};

struct point_3d {
    real_cpu x, y, z;
};

struct fiber_coords {
    real_cpu f[3];
    real_cpu s[3];
    real_cpu n[3];
};

struct fiber_coords_scale {
    real_cpu f[3];
    real_cpu s[3];
    real_cpu n[3];
    real_cpu x[3];
};

struct condutivity {
    real_cpu x, y, z, xy, xz, yz;
    struct fiber_coords fibers;
};

#define TRANSLATE(point, vx, vy, vz) POINT3D((point).x + (vx), (point).y + (vy), (point).z + (vz))
#define POINT3D(x, y, z)                                                                                                                                       \
    (struct point_3d) {                                                                                                                                        \
        x, y, z                                                                                                                                                \
    }
#define SAME_POINT3D(x)                                                                                                                                        \
    (struct point_3d) {                                                                                                                                        \
        x, x, x                                                                                                                                                \
    }

struct point_hash_entry {
    struct point_3d key;
    float value;
};

struct cell_hash_entry {
    struct point_3d key;
    struct point_3d value;
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

typedef uint32_t *ui32_array;
typedef struct element *element_array;
typedef float *f32_array;
typedef real *real_array;
typedef struct point_3d *point3d_array;
typedef struct line *line_array;
typedef int64_t *int64_array;
typedef uint8_t *ui8_array;
typedef int *int_array;
typedef char **string_array;

struct simulation_files {
    sds base_dir;
    string_array files_list;
    f32_array timesteps;
};

#define STRING_HASH_PRINT_KEY_VALUE(d)                                                                                                                         \
    do {                                                                                                                                                       \
        for(int64_t i__ = 0; i__ < shlen(d); i__++) {                                                                                                          \
            struct string_hash_entry e = (d[i__];                                                                                                              \
            printf("%s = %s\n", e.key, e.value);                                                                                                               \
        }                                                                                                                                                      \
    } while(0)

#define STRING_HASH_PRINT_KEY_VALUE_LOG(tag, d)                                                                                                                \
    do {                                                                                                                                                       \
        for(int64_t __i = 0; __i < shlen(d); __i++) {                                                                                                          \
            struct string_hash_entry __e = (d)[__i];                                                                                                           \
            log_info("%s %s = %s\n", tag, __e.key, __e.value);                                                                                                 \
        }                                                                                                                                                      \
    } while(0)

#define STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE(d, fn)                                                                                                 \
    do {                                                                                                                                                       \
        for(int64_t __i = 0; __i < hmlen(d); __i++) {                                                                                                          \
            struct string_voidp_hash_entry __e = (d)[__i];                                                                                                     \
            fn((struct config *)__e.value);                                                                                                                    \
        }                                                                                                                                                      \
    } while(0)

#define STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE_AND_KEY(d, fn)                                                                                         \
    do {                                                                                                                                                       \
        for(int64_t i__ = 0; i__ < hmlen(d); i__++) {                                                                                                          \
            struct string_voidp_hash_entry e = (d)[i__];                                                                                                       \
            fn(e.value, e.key);                                                                                                                                \
        }                                                                                                                                                      \
    } while(0)

#define STIM_CONFIG_HASH_FOR_INIT_FUNCTIONS(d)                                                                                                                 \
    do {                                                                                                                                                       \
        for(int64_t i__ = 0; i__ < hmlen(d); i__++) {                                                                                                          \
            struct string_voidp_hash_entry e = (d)[i__];                                                                                                       \
            init_config_functions((struct config *)e.value, "./shared_libs/libdefault_stimuli.so", e.key);                                                     \
        }                                                                                                                                                      \
    } while(0)

#define MODIFY_DOMAIN_CONFIG_HASH_FOR_INIT_FUNCTIONS(d)                                                                                                        \
    do {                                                                                                                                                       \
        for(int64_t __i = 0; __i < hmlen(d); __i++) {                                                                                                          \
            struct string_voidp_hash_entry e__ = (d)[__i];                                                                                                     \
            init_config_functions((struct config *)e__.value, "./shared_libs/libdefault_modify_domain.so", e__.key);                                           \
        }                                                                                                                                                      \
    } while(0)

#endif // MONOALG3D_COMMON_TYPES_H
