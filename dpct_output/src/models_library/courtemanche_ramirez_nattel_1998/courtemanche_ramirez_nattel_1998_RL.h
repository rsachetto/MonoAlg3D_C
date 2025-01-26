#ifndef MONOALG3D_MODEL_COURTEMANCHE_RAMIREZ_NATTEL_1998_RL_H
#define MONOALG3D_MODEL_COURTEMANCHE_RAMIREZ_NATTEL_1998_RL_H

#include "../model_common.h"

#define NEQ 21
#define INITIAL_V (-81.180000f)

#ifdef __CUDACC__

#include "../gpu_utils/gpu_utils.h"

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt, size_t pitch, bool use_adpt_dt);

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt, real min_dt);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, size_t pitch, bool use_adpt,
                          real abstol, real reltol, real max_dt);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt);
inline void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, struct ode_solver *solver);

void solve_model_ode_cpu(real dt, real *sv, real stim_current);

#endif //MONOALG3D_MODEL_COURTEMANCHE_RAMIREZ_NATTEL_1998_RL_H

