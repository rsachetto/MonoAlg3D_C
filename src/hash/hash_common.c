#include "hash_common.h"

bool point_equals(struct point_3d a, struct point_3d b) {
    return (a.x == b.x) && (a.y == b.y) && (a.z == b.z);
}


