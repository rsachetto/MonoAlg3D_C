#ifndef MONOALG3D_MIXED_FHN_MITCHELL_H
#define MONOALG3D_MIXED_FHN_MITCHELL_H

#include <stdint.h>
#include "model_common.h"

// Model 1 = Modified FitzHugh-Nagumo
// Model 2 = Mitchell-Shaeffer

#define NEQ_1 2         // Number of equations from model 1
#define NEQ_2 2         // Number of equations from model 2

#define INITIAL_V_1 (0.0)                       // Initial transmembrane potential from model 1
#define INITIAL_V_2 (0.00000820413566106744)    // Initial transmembrane potential from model 2

#ifdef __CUDACC__

extern "C" {
    #include "../utils/file_utils.h"
}

__constant__  size_t pitch;
size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(real *sv, uint32_t *mapping, int num_volumes);

__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t *mapping, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu_fhn(real *sv_, real *rDY_, real stim_current, int threadID_);
inline __device__ void RHS_gpu_mitchell(real *sv_, real *rDY_, real stim_current, int threadID_);

#else
#include "../utils/file_utils.h"
#endif

void solve_model_ode_cpu_fhn(real dt, real *sv, real stim_current);
void RHS_cpu_fhn(const real *sv, real *rDY_, real stim_current);

void solve_model_ode_cpu_mitchell(real dt, real *sv, real stim_current);
void RHS_cpu_mitchell(const real *sv, real *rDY_, real stim_current);

#endif //MONOALG3D_FHN_MOD_H

