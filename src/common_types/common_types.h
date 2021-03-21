#ifndef MONOALG3D_COMMON_TYPES_H
#define MONOALG3D_COMMON_TYPES_H

#include <stdint.h>

#define Pragma(x) _Pragma(#x)
#define OMP(directive) Pragma(omp directive)

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

#define REAL_DOUBLE

// Precision to be used for the calculations
#ifdef REAL_DOUBLE
typedef double real;
#else
typedef float real;
#endif

#define MALLOC_BYTES(type, bytes) (type *) malloc(bytes)
#define MALLOC_ONE_TYPE(type) (type *) malloc(sizeof(type))
#define MALLOC_ARRAY_OF_TYPE(type, n) (type *) malloc(sizeof(type) * n)

#define CALLOC_ONE_TYPE(type) (type *) calloc(1, sizeof(type))
#define CALLOC_ARRAY_OF_TYPE(type, n) (type *) calloc(n, sizeof(type))

typedef double real_cpu;

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

struct condutivity {
    real_cpu x, y, z, xy, xz, yz;
    struct fiber_coords fibers;
};

#define TRANSLATE(point, vx, vy, vz) POINT3D(point.x + vx, point.y + vy, point.z + vz)
#define POINT3D(x, y, z)                                                                                                                                       \
    (struct point_3d) {                                                                                                                                        \
        x, y, z                                                                                                                                                \
    }
#define SAME_POINT3D(x)                                                                                                                                        \
    (struct point_3d) {                                                                                                                                        \
        x, x, x                                                                                                                                                \
    }
#define ZERO_POINT3D POINT3D(0, 0, 0)

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
    string_array files_list;
    f32_array timesteps;
};

#define STRING_HASH_PRINT_KEY_VALUE(d)                                                                                                                         \
    do {                                                                                                                                                       \
        for(long i = 0; i < shlen(d); i++) {                                                                                                                   \
            struct string_hash_entry e = d[i];                                                                                                                 \
            printf("%s = %s\n", e.key, e.value);                                                                                                               \
        }                                                                                                                                                      \
    } while(0)

#define STRING_HASH_PRINT_KEY_VALUE_LOG(d)                                                                                                                     \
    do {                                                                                                                                                       \
        for(long i = 0; i < shlen(d); i++) {                                                                                                                   \
            struct string_hash_entry e = d[i];                                                                                                                 \
            log_info("%s = %s\n", e.key, e.value);                                                                                               \
        }                                                                                                                                                      \
    } while(0)

#define STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE(d, fn)                                                                                                 \
    do {                                                                                                                                                       \
        for(long i = 0; i < hmlen(d); i++) {                                                                                                                   \
            struct string_voidp_hash_entry e = d[i];                                                                                                           \
            fn(e.value);                                                                                                                                       \
        }                                                                                                                                                      \
    } while(0)

#define STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE_AND_KEY(d, fn)                                                                                         \
    do {                                                                                                                                                       \
        for(long i = 0; i < hmlen(d); i++) {                                                                                                                   \
            struct string_voidp_hash_entry e = d[i];                                                                                                           \
            fn(e.value, e.key);                                                                                                                                \
        }                                                                                                                                                      \
    } while(0)

#define STIM_CONFIG_HASH_FOR_INIT_FUNCTIONS(d)                                                                                                                 \
    do {                                                                                                                                                       \
        for(long i = 0; i < hmlen(d); i++) {                                                                                                                   \
            struct string_voidp_hash_entry e = d[i];                                                                                                           \
            init_config_functions(e.value, "./shared_libs/libdefault_stimuli.so", e.key);                                                                      \
        }                                                                                                                                                      \
    } while(0)

#define MODIFY_DOMAIN_CONFIG_HASH_FOR_INIT_FUNCTIONS(d)                                                                                                        \
    do {                                                                                                                                                       \
        for(long i = 0; i < hmlen(d); i++) {                                                                                                                   \
            struct string_voidp_hash_entry e = d[i];                                                                                                           \
            init_config_functions(e.value, "./shared_libs/libdefault_modify_domain.so", e.key);                                                                \
        }                                                                                                                                                      \
    } while(0)

#endif
