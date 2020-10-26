#ifndef MONOALG3D_MODEL_ARPF_2009_H
#define MONOALG3D_MODEL_ARPF_2009_H

#include "model_common.h"

// MODEL INFO: Mathematical models of the electrical action potential of Purkinje fibre cells, 2009

#define NEQ 20
#define INITIAL_V (-69.1370441635924)

#ifdef __CUDACC__

#include "../gpu_utils/gpu_utils.h"

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt);
inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id);

#endif

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt);
inline void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int thread_id);


#endif //MONOALG3D_FHN_MOD_H

