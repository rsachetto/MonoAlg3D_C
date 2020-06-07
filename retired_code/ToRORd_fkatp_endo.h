#ifndef MONOALG3D_MODEL_TORORD_FKATP_ENDO_H
#define MONOALG3D_MODEL_TORORD_FKATP_ENDO_H

#include "model_common.h"

#define NEQ 43
#define INITIAL_V (-88.763800f)

#ifdef __CUDACC__

__constant__  size_t pitch;
size_t pitch_h;

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes);
__global__ void kernel_set_model_initial_conditions_matlab(real *sv, int num_volumes);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt);
inline __device__ void RHS_gpu_matlab(real *sv_, real *rDY_, real stim_current, int threadID_, real dt);
__inline__ __device__ void getINa_Grandi(real *INa, real *dm, real *dh, real *dhp, real *dj, real *djp, real v, real m, real h, real hp, real j, real jp, real fINap, real ENa, real INa_Multiplier);
__inline__ __device__  void getINaL_ORd2011(real *INaL, real *dmL, real *dhL, real *dhLp, real v, real mL, real hL, real hLp, real fINaLp, real ENa, int celltype, real INaL_Multiplier);
__inline__ __device__ void getICaL_ORd2011_jt(real *ICaL_ss, real *ICaNa_ss, real *ICaK_ss, real *ICaL_i, real *ICaNa_i, real *ICaK_i,
                                              real *dd, real *dff, real *dfs, real *dfcaf,real * dfcas,real * djca, real *dnca, real * dnca_i,
                                              real *dffp, real *dfcafp, real *PhiCaL_ss, real *PhiCaL_i, real *gammaCaoMyo, real *gammaCaiMyo,
                                              real v, real d, real ff, real fs, real fcaf, real fcas, real jca, real nca, real nca_i, real ffp,
                                              real fcafp, real fICaLp, real cai, real cass, real cao, real nai, real nass, real nao, real ki,
                                              real kss, real ko, real cli, real clo, real celltype, real ICaL_fractionSS, real ICaL_PCaMultiplier);

__inline__ __device__ void getIKr_ORd2011_MM(real *IKr, real *dc0, real *dc1, real *dc2, real *do_, real *di, real V, real c0,
                                             real c1, real c2, real  o, real i, real ko, real EK, int celltype, real IKr_Multiplier);

__inline__ __device__ void getIKs_ORd2011(real *IKs, real *dxs1, real *dxs2, real v, real xs1, real xs2, real cai, real EKs, int celltype, real IKs_Multiplier);

__inline__ __device__ void getIK1_CRLP(real *IK1, real v, real ko, real  EK, int celltype, real IK1_Multiplier);

__inline__ __device__ void getINaCa_ORd2011(real *INaCa_i, real *INaCa_ss, real v, real F, real R, real T, real nass,
                                            real nai, real nao, real cass, real cai, real cao, int celltype,
                                            real INaCa_Multiplier, real INaCa_fractionSS);

__inline__ __device__ void getINaK_ORd2011(real *INaK, real v, real F, real R, real T, real nai, real nao, real ki, real ko, int celltype, real  INaK_Multiplier);

__inline__ __device__ void getJrel_ORd2011(real *Jrel, real *dJrelnp, real *dJrelp, real Jrelnp, real Jrelp, real ICaL, real cass, real cajsr, real fJrelp, int celltype, real Jrel_Multiplier);

__inline__ __device__ void getJup_ORd2011(real *Jup, real *Jleak, real cai, real cansr, real fJupp, int celltype, real Jup_Multiplier);

__inline__ __device__ void getITo_ORd2011(real *Ito, real *da, real *diF, real *diS, real *dap, real *diFp, real *diSp, real v, real a,
                                          real iF, real iS, real ap, real iFp, real iSp, real fItop, real EK, int celltype, real Ito_Multiplier);


#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt);

void solve_model_ode_cpu(real dt, real *sv, real stim_current, real **Z, real **F, int *count, int *m, int sv_id);
#endif //MONOALG3D_MODEL_TORORD_FKATP_ENDO_H

