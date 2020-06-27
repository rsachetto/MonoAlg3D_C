#ifndef MONOALG3D_MIXED_TT2004_MYO_EPI_H
#define MONOALG3D_MIXED_TT2004_MYO_EPI_H

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "../../model_common.h"

// Model 1 = TenTusscher2004 MCELL
// Model 2 = TenTusscher2004 EPI

#define NEQ 17         

#define EPI

#define INITIAL_V (-86.2)                       // Initial transmembrane potential

#ifdef __CUDACC__

__constant__  size_t pitch;
size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(real *sv, real*IC, uint32_t *mapping, int num_volumes);

__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps,
                          uint32_t *mapping);

inline __device__ void RHS_gpu_myo (real *sv, real *rDY_, real stim_current, int threadID_, real dt);
inline __device__ void RHS_gpu_epi (real *sv, real *rDY_, real stim_current, int threadID_, real dt);

#endif

void solve_model_ode_cpu_myo (real dt, real *sv, real stim_current);
void RHS_cpu_myo (const real *sv, real *rDY_, real stim_current, real dt);

void solve_model_ode_cpu_epi (real dt, real *sv, real stim_current);
void RHS_cpu_epi (const real *sv, real *rDY_, real stim_current, real dt);

#endif //MONOALG3D_MIXED_TT2004_MYO_EPI_H

