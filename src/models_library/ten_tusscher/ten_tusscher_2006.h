#ifndef MONOALG3D_MODEL_TEN_TUSSCHER_2006_H
#define MONOALG3D_MODEL_TEN_TUSSCHER_2006_H

#include "../model_common.h"

#define NEQ 19
#define INITIAL_V (-85.23f)

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

__global__ void kernel_set_model_inital_conditions(real *sv, real *initial_conditions, int num_volumes);

__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt);
void solve_model_ode_cpu(real dt, real *sv, real stim_current);

#endif //MONOALG3D_MODEL_TEN_TUSSCHER_2006_H
