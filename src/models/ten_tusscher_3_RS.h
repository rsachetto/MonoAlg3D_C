#ifndef MONOALG3D_MODEL_TEN_TUSSCHER_COMMON_H
#define MONOALG3D_MODEL_TEN_TUSSCHER_COMMON_H

#include <unitypes.h>
#include "../main/constants.h"

#define ENDO
#define NEQ 12
#define INITAL_V -85.23f

#ifdef __CUDACC__
static __device__ size_t pitch;
static size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(Real *sv, int num_volumes);

__global__ void solve_gpu(Real dt, Real *sv, Real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, Real *fibrosis, Real atpi);

inline __device__ void RHS_gpu(Real *sv_, Real *rDY_, Real stim_current, int threadID_, Real dt, Real fibrosis, Real atpi);
#else
void RHS_cpu(const Real *sv, Real *rDY_, Real stim_current, Real dt, Real fibrosis, Real atpi);
void solve_model_ode_cpu(Real dt, Real *sv, Real stim_current, Real fibrosis, Real atpi );
#endif

#endif //MONOALG3D_MODEL_TEN_TUSSCHER_COMMON_H
