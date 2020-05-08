#ifndef MONOALG3D_MODEL_COURTEMANCHE_RAMIREZ_NATTEL_1998_H
#define MONOALG3D_MODEL_COURTEMANCHE_RAMIREZ_NATTEL_1998_H

#include "model_common.h"
#include "model_gpu_utils.h"

#define NEQ 21
#define INITIAL_V (-81.180000f)

#ifdef __CUDACC__

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int thread_id, real dt);
inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id);

#endif

void RHS_cpu(const real *sv, real *rDY, real stim_current, real dt);
inline void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int thread_id);

void solve_model_ode_cpu(real dt, real *sv, real stim_current);

#endif //MONOALG3D_MODEL_COURTEMANCHE_RAMIREZ_NATTEL_1998_H

