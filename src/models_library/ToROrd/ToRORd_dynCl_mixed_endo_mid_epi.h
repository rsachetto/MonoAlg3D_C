#ifndef MONOALG3D_MODEL_TORORD_DYNCL_MIXED_ENDO_MID_EPI_H
#define MONOALG3D_MODEL_TORORD_DYNCL_MIXED_ENDO_MID_EPI_H

// TOMEK, Jakub and BUENO-OROVIO, Alfonso and RODRIGUEZ, Blanca  
// ToR-ORd-dynCl: an update of the ToR-ORd model of human ventricular cardiomyocyte with dynamic intracellular chloride 
//  bioRxiv, 2020.

#include "../model_common.h"
#include "../../extra_data_library/helper_functions.h"

#define NEQ 45
#define INITIAL_V (-9.035192e+01)

#define ENDO 0.0
#define MID  1.0
#define EPI  2.0

// CPU macros
#define SOLVE_EQUATION_EULER_CPU(id) sv[id] = dt * rDY[id] + rY[id]

#define SOLVE_EQUATION_RUSH_LARSEN_CPU(id) sv[id] = (abs(a[id]) < TOLERANCE) ? rY[id] + dt * (rY[id] * a[id] + b[id]) : \
                                        exp(a[id] * dt)*(rY[id] + (b[id] / a[id])) - (b[id] / a[id] )

#define SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(id) edos_old_aux_[id] = sv[id]; \
                                                 edos_new_euler_[id] = _k1__[id] * *dt + edos_old_aux_[id]; \
                                                 sv[id] = edos_new_euler_[id]

#define SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(id) edos_old_aux_[id] = sv[id]; \
                                                 edos_new_euler_[id] = (fabs(a_[id]) < abs_tol) ? edos_old_aux_[id] + (edos_old_aux_[id] * a_[id] + b_[id])*(*dt) : exp(a_[id]*(*dt))*(edos_old_aux_[id] + (b_[id] / a_[id])) - (b_[id] / a_[id]); \
                                                 sv[id] = edos_new_euler_[id]

#define SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(id) _k2__[id] = rDY[id]; \
                                                    f = (_k1__[id] + _k2__[id]) * 0.5; \
                                                    y_2nd_order = edos_old_aux_[id] + (*dt) * f; \
                                                    auxError = (fabs(y_2nd_order) < abs_tol) ? fabs(edos_new_euler_[id] - abs_tol) : fabs( (y_2nd_order - edos_new_euler_[id])/(y_2nd_order) ); \
                                                    greatestError = (auxError > greatestError) ? auxError : greatestError

#define SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(id)    _k2__[id] = rDY[id]; \
                                                    as = (a_[id] + a_new[id]) * 0.5; \
                                                    bs = (b_[id] + b_new[id]) * 0.5; \
                                                    y_2nd_order = (fabs(as) < abs_tol) ? edos_old_aux_[id] + (*dt) * (edos_old_aux_[id]*as + bs) : exp(as*(*dt))*(edos_old_aux_[id] + (bs/as)) - (bs/as); \
                                                    auxError = (fabs(y_2nd_order) < abs_tol) ? fabs(edos_new_euler_[id] - abs_tol) : fabs( (y_2nd_order - edos_new_euler_[id])/(y_2nd_order) ); \
                                                    greatestError = (auxError > greatestError) ? auxError : greatestError

// GPU macros
#define SOLVE_EQUATION_EULER_GPU(id) *((real *)((char *)sv + pitch * id) + sv_id) = *((real *)((char *)sv + pitch * id) + sv_id) + dt * rDY[id]

#define SOLVE_EQUATION_RUSH_LARSEN_GPU(id) *((real *)((char *)sv + pitch * id) + sv_id) = (abs(a[id]) < TOLERANCE) ? *((real *)((char *)sv + pitch * id) + sv_id) + dt * ( *((real *)((char *)sv + pitch * id) + sv_id) * a[id] + b[id]) : \
                                                                                    exp(a[id] * dt)*(*((real *)((char *)sv + pitch * id) + sv_id) + (b[id] / a[id])) - (b[id] / a[id]) 

#define SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(id) edos_old_aux_[id] = sv_local[id]; \
                                                    edos_new_euler_[id] = _k1__[id] * dt + edos_old_aux_[id]; \
                                                    sv_local[id] = edos_new_euler_[id]

#define SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(id) edos_old_aux_[id] = sv_local[id]; \
                                                    edos_new_euler_[id] = (a_[id] < abstol) ? edos_old_aux_[id] + (edos_old_aux_[id] * a_[id] + b_[id])*(dt) : exp(a_[id]*(dt))*(edos_old_aux_[id] + (b_[id] / a_[id])) - (b_[id] / a_[id]); \
                                                    sv_local[id] = edos_new_euler_[id]

#define SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(id) _k2__[id] = rDY[id]; \
                                                    f = (_k1__[id] + _k2__[id]) * 0.5; \
                                                    y_2nd_order = edos_old_aux_[id] + (dt) * f; \
                                                    auxError = (fabs(y_2nd_order) < abstol) ? fabs(edos_new_euler_[id] - abstol) : fabs( (y_2nd_order - edos_new_euler_[id])/(y_2nd_order) ); \
                                                    greatestError = (auxError > greatestError) ? auxError : greatestError

#define SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(id)    _k2__[id] = rDY[id]; \
                                                    as = (a_[id] + a_new[id]) * 0.5; \
                                                    bs = (b_[id] + b_new[id]) * 0.5; \
                                                    y_2nd_order = (fabs(as) < abstol) ? edos_old_aux_[id] + (dt) * (edos_old_aux_[id]*as + bs) : exp(as*(dt))*(edos_old_aux_[id] + (bs/as)) - (bs/as); \
                                                    auxError = (fabs(y_2nd_order) < abstol) ? fabs(edos_new_euler_[id] - abstol) : fabs( (y_2nd_order - edos_new_euler_[id])/(y_2nd_order) ); \
                                                    greatestError = (auxError > greatestError) ? auxError : greatestError

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt);
__global__ void kernel_set_model_initial_conditions_endo_mid_epi(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt,\
                                                real *initial_endo, real *initial_epi, real *initial_mid, real *transmurality);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, real *extra_params,\
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt);
__global__ void solve_endo_mid_epi_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, real *transmurality, real *extra_params,\
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt);

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, real transmurality, real *extra_params, int threadID_, real dt, size_t pitch, bool use_adpt_dt);
inline __device__ void RHS_RL_gpu(real *a_, real *b_, real *sv, real *rDY_, real stim_current, real transmurality, real *extra_params, int threadID_, real dt, size_t pitch, bool use_adpt_dt);
inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real transmurality, real *extra_params, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt);
inline __device__ void solve_rush_larsen_gpu_adpt(real *sv, real stim_curr, real transmurality, real *extra_params, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real const *extra_params);
void RHS_RL_cpu (real *a_, real *b_, const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real const *extra_params);
void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real transmurality, real final_time, int sv_id, struct ode_solver *solver, real const *extra_params);
void solve_rush_larsen_cpu_adpt(real *sv, real stim_curr, real transmurality, real final_time, int sv_id, struct ode_solver *solver, real const *extra_params);
void solve_model_ode_cpu(real dt, real *sv, real stim_current, real transmurality, real const *extra_params);

#endif //MONOALG3D_MODEL_TORORD_DYNCL_MIXED_ENDO_MID_EPI_H

