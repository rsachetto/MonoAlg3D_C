#ifndef MONOALG3D_MODEL_MINIMAL_MODEL_H
#define MONOALG3D_MODEL_MINIMAL_MODEL_H

#include "../model_common.h"
#include "../../extra_data_library/helper_functions.h"

/*
Minimal model for human ventricular action potentials in tissue,
Alfonso Bueno-Orovio, Elizabeth M Cherry, Flavio H Fenton,
(2008)
*/

#define NEQ 4
#define INITIAL_V (0.0f)

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

static __device__ size_t pitch;
static size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(real *sv, real* ICs, int num_volumes);

__global__ void solve_gpu(real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps,uint32_t *mapping);

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt,
                                int type_cell);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, int type_cell);
void solve_model_ode_cpu(real dt, real *sv, real stim_current, int type_cell);

#endif