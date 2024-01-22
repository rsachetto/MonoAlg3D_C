#ifndef MONOALG3D_MODEL_CRN2017_TASK1_IKACh_H
#define MONOALG3D_MODEL_CRN2017_TASK1_IKACh_H

#include "../model_common.h"

#define NEQ 26
#define INITIAL_V (-82.7f)

#ifdef __CUDACC__

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, real *extra_parameters);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, real *extra_parameters);

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int threadID, real dt, real *extra_parameters, int mapping);
inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int threadID, real *extra_parameters, int mapping);

#endif

void RHS_cpu(const real *sv, real *rDY, real stim_current, real dt, real *extra_parameters, int mapping);
inline void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int threadID, real *extra_parameters, int mapping);

void solve_model_ode_cpu(real dt, real *sv, real stim_current, real *extra_parameters, int mapping);

#endif //MONOALG3D_MODEL_CRN2017_TASK1_IKACh_H

