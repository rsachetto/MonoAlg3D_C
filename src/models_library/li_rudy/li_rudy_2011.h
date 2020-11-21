#ifndef MONOALG3D_MODEL_LIRUDY_2011_H
#define MONOALG3D_MODEL_LIRUDY_2011_H

#include "../model_common.h"
#include <stdint.h>

#define NEQ 41
#define INITIAL_V (-84.058830)

#ifdef __CUDACC__

__constant__  size_t pitch;
size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes);

__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, real dt, int threadID_);

#endif

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(const real *sv, real *rDY_, real dt, real stim_current);

#endif // MONOALG3D_MODEL_LIRUDY_2011_H

