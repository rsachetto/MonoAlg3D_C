//
// Created by sachetto on 29/04/2020.
//

#ifndef MONOALG3D_C_MODEL_COMMON_H
#define MONOALG3D_C_MODEL_COMMON_H

#include "../common_types/common_types.h"
#include "../ode_solver/ode_solver.h"
#include <float.h>

#ifdef CELL_MODEL_REAL_DOUBLE
    #ifdef __CUDACC__
        #if __CUDACC_VER_MAJOR__ > 9
            #pragma message "Using double precision for the GPU cellular model"
            #include "set_double_precision.h"
        #else
            #pragma message "Using single precision for the GPU cellular model"
            #include "set_single_precision.h"
        #endif
    #else
    #pragma message "Using double precision for the CPU cellular model"
        #include "set_double_precision.h"
    #endif
#else
    #ifdef __CUDACC__
        #pragma message "Using single precision for the GPU cellular model"
    #else
        #pragma message "Using single precision for the CPU cellular model"
    #endif
    #include "set_single_precision.h"
#endif

#endif // MONOALG3D_C_MODEL_COMMON_H
