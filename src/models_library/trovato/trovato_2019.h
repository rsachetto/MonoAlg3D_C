#ifndef MONOALG3D_MODEL_TROVATO_2019_H
#define MONOALG3D_MODEL_TROVATO_2019_H

#include "../model_common.h"

#define NEQ 46
#define INITIAL_V (-86.7099)

#include "../default_solvers.h"

#ifdef __CUDACC__

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt);
inline __device__ void solve_rush_larsen_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt);
inline __device__ void RHS_RL_gpu(real *a_, real *b_, real *sv, real *rDY_, real stim_current, int threadID_, real dt, size_t pitch, bool use_adpt_dt);

#endif

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver);
void solve_rush_larsen_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver);
void RHS_RL_cpu (real *a_, real *b_, real *sv, real *rDY_, real stim_current, real dt);

#endif //MONOALG3D_MODEL_TROVATO_2019_H

