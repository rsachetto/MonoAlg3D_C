//
// Created by sachetto on 29/04/2020.
//

#ifndef MONOALG3D_C_MODEL_COMMON_H
#define MONOALG3D_C_MODEL_COMMON_H

#include "../common_types/common_types.h"
#include "../ode_solver/ode_solver.h"
#include <float.h>

#ifndef GPU_REAL_DOUBLE
#pragma message "Using single precision for the GPU model"
#define pow powf
#define sqrt sqrtf
#define exp expf
#define fabs fabsf
#define log logf
#define REAL_MIN FLT_MIN
#define REAL_MAX FLT_MAX
#define REAL_EPS FLT_EPSILON
#else
#pragma message "Using double precision for the GPU model"
#define REAL_MIN DBL_MIN
#define REAL_MAX DBL_MAX
#define REAL_EPS DBL_EPSILON
#endif

#endif // MONOALG3D_C_MODEL_COMMON_H
