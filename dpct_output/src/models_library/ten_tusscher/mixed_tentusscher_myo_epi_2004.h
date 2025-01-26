#ifndef MONOALG3D_MIXED_MODEL_TEN_TUSSCHER_MYO_EPI_2004_H
#define MONOALG3D_MIXED_MODEL_TEN_TUSSCHER_MYO_EPI_2004_H

// Model 1 = TenTusscher 2004 myocardium
// Model 2 = TenTusscher 2004 epicardium

#include "../model_common.h"

#define NEQ 17
#define INITIAL_V (-86.2f)

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

extern "C" {
    #include "../../utils/file_utils.h"
}

__global__ void kernel_set_model_inital_conditions(real *sv, uint32_t *mapping, int num_volumes);

__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t *mapping, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu_myo(real *sv_, real *rDY_, real stim_current, int threadID_, real dt);
inline __device__ void RHS_gpu_epi(real *sv_, real *rDY_, real stim_current, int threadID_, real dt);

#else
#include "../../utils/file_utils.h"
#endif

void solve_model_ode_cpu_myo(real dt, real *sv, real stim_current);
void RHS_cpu_myo(const real *sv, real *rDY_, real stim_current, real dt);

void solve_model_ode_cpu_epi(real dt, real *sv, real stim_current);
void RHS_cpu_epi(const real *sv, real *rDY_, real stim_current, real dt);

#endif // MONOALG3D_MIXED_MODEL_TEN_TUSSCHER_MYO_EPI_2004_H
