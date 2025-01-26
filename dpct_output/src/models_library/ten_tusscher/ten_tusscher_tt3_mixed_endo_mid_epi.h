#ifndef MONOALG3D_MODEL_TEN_TUSSCHER_3_MIXED_ENDO_MID_EPI_H
#define MONOALG3D_MODEL_TEN_TUSSCHER_3_MIXED_ENDO_MID_EPI_H

#include "../model_common.h"
#include "../../extra_data_library/helper_functions.h"

#define NEQ 12
#define INITIAL_V (-86.2f)

#define ENDO 0.0
#define MID  1.0
#define EPI  2.0

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

static __device__ size_t pitch;
static size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(real *sv, real* ICs, int num_volumes);

__global__ void solve_gpu(real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, real *fibrosis, real *transmurality, real *extra_parameters);

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt, real fibrosis, real transmurality,
                               real *extra_parameters);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real fibrosis, real transmurality, real const *extra_parameters);
void solve_model_ode_cpu(real dt, real *sv, real stim_current, real fibrosis, real transmurality, real *extra_parameters);

#endif // MONOALG3D_MODEL_TEN_TUSSCHER_3_MIXED_ENDO_MID_EPI_H
