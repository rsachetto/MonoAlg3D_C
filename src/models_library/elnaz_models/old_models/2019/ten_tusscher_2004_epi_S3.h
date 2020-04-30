//AP+Rd+Vmax:3obj: Tr2,pop41,
//Error:0.237161115433174,0.000785891679072392,0.000207784955952597 (calculated Rc=0.00523558046455782)
//GNa,GNab,GCaL,GCab,Gto,GKr,GKs,GK1,Gpk,PNak,KNaCa,Vmax_up,GpCa,arel,crel,Vleak,
//parameters: 14.4350685070016	3.57702620495540e-05	0.000142934187751013	0.000433192953639154	0.305652893755267	0.141573793892399	0.209159988702282	4.93898438176732	0.0172111589141696	1.54854638201186	1099.86953805106	0.000467408934683950	0.217083954885964	0.00987840511952626	0.00571383652914274	3.60633195830624e-05

#ifndef MONOALG3D_MODEL_TEN_TUSSCHER_2004_H
#define MONOALG3D_MODEL_TEN_TUSSCHER_2004_H

#include <stdint.h>
#include "model_common.h"

#define NEQ 17
#define INITIAL_V (-86.2f)
#define EPI

#ifdef __CUDACC__

extern "C" {
    #include "../utils/file_utils.h"
}

__constant__  size_t pitch;
size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes);

__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt);

#else
#include "../utils/file_utils.h"
#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt);
void solve_model_ode_cpu(real dt, real *sv, real stim_current);

#endif //MONOALG3D_MODEL_TEN_TUSSCHER_2004_H
