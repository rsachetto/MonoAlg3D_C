#ifndef MONOALG3D_MODEL_TORORD_FKATP_MIXED_ENDO_MID_EPI_H
#define MONOALG3D_MODEL_TORORD_FKATP_MIXED_ENDO_MID_EPI_H

// TOMEK, Jakub et al. Development, calibration, and validation of a novel human ventricular myocyte model in health, disease, and drug block. 
//  Elife, v. 8, p. e48890, 2019.

#include "../model_common.h"

#define NEQ 43
#define INITIAL_V (-88.763800)

#define ENDO 0.0
#define MID  1.0
#define EPI  2.0

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt);
__global__ void kernel_set_model_initial_conditions_endo_mid_epi(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt,\
                                                real *initial_endo, real *initial_epi, real *initial_mid, real *mapping);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve,\
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt);
__global__ void solve_endo_mid_epi_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, real *mapping,\
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt);

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, real mapping, int threadID_, real dt, size_t pitch, bool use_adpt_dt);
inline __device__ void RHS_RL_gpu(real *a_, real *b_, real *sv, real *rDY_, real stim_current, real mapping, int threadID_, real dt, size_t pitch, bool use_adpt_dt);
inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt);
inline __device__ void solve_rush_larsen_gpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real mapping);
void RHS_RL_cpu (real *a_, real *b_, const real *sv, real *rDY_, real stim_current, real dt, real mapping);
void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int sv_id, struct ode_solver *solver);
void solve_rush_larsen_cpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int sv_id, struct ode_solver *solver);
void solve_model_ode_cpu(real dt, real *sv, real stim_current, real mapping);

#endif //MONOALG3D_MODEL_TORORD_FKATP_MIXED_ENDO_MID_EPI_H

