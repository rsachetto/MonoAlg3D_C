//
// Created by sachetto on 01/10/17.
//

#ifndef MONOALG3D_UTILS_H
#define MONOALG3D_UTILS_H

#include <stddef.h>

#include "../common_types/common_types.h"

void sort_vector(real_cpu **a, int length);
int inside_mesh(real_cpu **a, real_cpu x, real_cpu y, real_cpu z, int first, int last);
float calculate_mean (float *arr, unsigned long size);

#endif //MONOALG3D_UTILS_H_H
