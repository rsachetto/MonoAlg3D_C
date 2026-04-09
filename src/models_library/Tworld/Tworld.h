#ifndef MONOALG3D_MODEL_TWORLD_H
#define MONOALG3D_MODEL_TWORLD_H

/*
    "T-World Virtual Human Cardiomyocyte. II. Organ-Scale Simulations and Applications"
    - Original implementation: 02/11/2024, Lucas Arantes Berg, University of Oxford
    - Reference: "T-World Virtual Human Cardiomyocyte. II. Organ-Scale Simulations and Applications"
        Journal = Circulation Research
	Year = 2026
	doi: https://doi.org/10.1161/CIRCRESAHA.125.328123
*/

#include "../model_common.h"

#define NEQ 93
#define INITIAL_V (-86.98548236471670236369)

#define ENDO 0
#define EPI  1
#define MID  2

// CPU macros
#define SOLVE_EQUATION_EULER_CPU(id) sv[id] = dt * rDY[id] + rY[id]

#define SOLVE_EQUATION_RUSH_LARSEN_CPU(id) sv[id] = (abs(a[id]) < TOLERANCE) ? rY[id] + dt * (rY[id] * a[id] + b[id]) : \
                                        exp(a[id] * dt)*(rY[id] + (b[id] / a[id])) - (b[id] / a[id] )

// GPU macros
#define SOLVE_EQUATION_EULER_GPU(id) *((real *)((char *)sv + pitch * id) + sv_id) = *((real *)((char *)sv + pitch * id) + sv_id) + dt * rDY[id]

#define SOLVE_EQUATION_RUSH_LARSEN_GPU(id) *((real *)((char *)sv + pitch * id) + sv_id) = (abs(a[id]) < TOLERANCE) ? *((real *)((char *)sv + pitch * id) + sv_id) + dt * ( *((real *)((char *)sv + pitch * id) + sv_id) * a[id] + b[id]) : \
                                                                                    exp(a[id] * dt)*(*((real *)((char *)sv + pitch * id) + sv_id) + (b[id] / a[id])) - (b[id] / a[id]) 

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt, size_t pitch, bool use_adpt_dt,
													 int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity);
													 
inline __device__ void RHS_RL_gpu (real *a_, real *b_, real *sv, real *rDY, real stim_current, int thread_id, real dt, size_t pitch, bool use_adpt, 
													 int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity);

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt, real min_dt, real *extra_data);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, size_t pitch, bool use_adpt,
                          real abstol, real reltol, real max_dt, real *extra_data);

#endif

void RHS_cpu (const real *sv, real *rDY, real stim_current, real dt, int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity);
void RHS_RL_cpu (real *a_, real *b_, const real *sv, real *rDY, real stim_current, real dt, int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity);
inline void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, struct ode_solver *solver, int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity);

void solve_model_ode_cpu(real dt, real *sv, real stim_current, int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity);

#endif //MONOALG3D_MODEL_TWORLD_H
