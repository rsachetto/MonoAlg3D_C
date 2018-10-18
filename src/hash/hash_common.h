#ifndef MONOALG3D_HASH_COMMON_H
#define MONOALG3D_HASH_COMMON_H

#define INITIAL_SIZE (8)
#define GROWTH_FACTOR (2)
#define MAX_LOAD_FACTOR (1)

#include <stdbool.h>

struct point_3d {
    float x, y, z;
};

bool point_equals(struct point_3d a, struct point_3d b);

#endif