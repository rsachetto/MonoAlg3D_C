#ifndef MONOALG3D_MODEL_TORORD_DYNCL_MIXED_ENDO_MID_EPI_H
#define MONOALG3D_MODEL_TORORD_DYNCL_MIXED_ENDO_MID_EPI_H

// TOMEK, Jakub and BUENO-OROVIO, Alfonso and RODRIGUEZ, Blanca  
// ToR-ORd-dynCl: an update of the ToR-ORd model of human ventricular cardiomyocyte with dynamic intracellular chloride 
//  bioRxiv, 2020.

#include "../model_common.h"

#define NEQ 45
#define INITIAL_V (-8.974808e+01)

#define ENDO 0.0
#define MID  1.0
#define EPI  2.0

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

__global__ void kernel_set_model_initial_conditions(real *sv,\
                                                real *initial_endo, real *initial_epi, real *initial_mid,\
                                                real *mapping, int num_volumes);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, real *mapping);

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, real mapping, int thread_id, real dt);
inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int thread_id);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real mapping);
inline void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int thread_id);

void solve_model_ode_cpu(real dt, real *sv, real stim_current, real mapping);

#endif //MONOALG3D_MODEL_TORORD_DYNCL_MIXED_ENDO_MID_EPI_H

