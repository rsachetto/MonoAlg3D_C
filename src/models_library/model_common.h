//
// Created by sachetto on 05/10/17.
//

// Every model need to implement the functions described in this model file in order to be loaded correctly from the
// edo solver. This models_library should compile without using any dependency of our codebase

#ifndef MONOALG3D_MODEL_COMMON_H
#define MONOALG3D_MODEL_COMMON_H

#include "../monodomain/constants.h"

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

struct cell_model_data {
    int number_of_ode_equations;
    real initial_v;
    char *model_library_path;
};

#define GET_CELL_MODEL_DATA(name) EXPORT_FN void name (struct cell_model_data *cell_model, bool get_initial_v, bool get_neq)
typedef GET_CELL_MODEL_DATA (get_cell_model_data_fn);

// CPU FUNCTIONS
#define SET_ODE_INITIAL_CONDITIONS_CPU(name) EXPORT_FN void name (real *sv)
typedef SET_ODE_INITIAL_CONDITIONS_CPU (set_ode_initial_conditions_cpu_fn);

#define SOLVE_MODEL_ODES_CPU(name)                                                                                     \
EXPORT_FN void name (real dt, real *sv, real *stim_currents, const uint32_t *cells_to_solve, uint32_t num_cells_to_solve,          \
               int num_steps, void *extra_data)
typedef SOLVE_MODEL_ODES_CPU (solve_model_ode_cpu_fn);

// GPU FUNCTIONS
#define SET_ODE_INITIAL_CONDITIONS_GPU(name) EXPORT_FN size_t name (real **sv, uint32_t num_volumes)
typedef SET_ODE_INITIAL_CONDITIONS_GPU (set_ode_initial_conditions_gpu_fn);

#define SOLVE_MODEL_ODES_GPU(name)                                                                                     \
EXPORT_FN void name (real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, uint32_t num_cells_to_solve,          \
               int num_steps, void *extra_data, size_t extra_data_bytes_size)
typedef SOLVE_MODEL_ODES_GPU(solve_model_ode_gpu_fn);

// typedef void (*update_gpu_fn_pt)(real *, uint32_t *, size_t , int );

////TODO: DEBUG: REMOVE
//__global__ void print_svs(real *sv, int numCells) {
//
//    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
//    if(threadID == 0) {
//        for (int i = 0; i < numCells; ++i) {
//            for (int j = 0; j < NEQ; ++j) {
//                printf("%lf\n", *((real*)((char*)sv + pitch * j)+i));
//
//            }
//
//        }
//    }
//}

#endif // MONOALG3D_MODEL_COMMON_H
