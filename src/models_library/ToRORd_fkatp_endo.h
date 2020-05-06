#ifndef MONOALG3D_MODEL_TORORD_FKATP_ENDO_H
#define MONOALG3D_MODEL_TORORD_FKATP_ENDO_H

#include "model_common.h"
#include "model_gpu_utils.h"

#define NEQ 43
#define INITIAL_V (-88.763800f)

#define IFNUMBER_1(name)if((celltype==1.000000000000000e+00f)) { (name) = (cmdnmax_b*1.300000000000000e+00f);    }  else{ (name) = cmdnmax_b;    }
#define IFNUMBER_2(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = 0.000000000000000e+00f;    }  else{ (name) = (5.700000000000000e-02f*expf(((-(v_old_+8.000000000000000e+01f))/6.800000000000000e+00f)));    }
#define IFNUMBER_3(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = (7.700000000000000e-01f/(1.300000000000000e-01f*(1.000000000000000e+00f+expf(((-(v_old_+1.066000000000000e+01f))/1.110000000000000e+01f)))));    }  else{ (name) = ((2.700000000000000e+00f*expf((7.900000000000000e-02f*v_old_)))+(3.100000e+05*expf((3.485000000000000e-01f*v_old_))));    }
#define IFNUMBER_4(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = 0.000000000000000e+00f;    }  else{ (name) = (((((-2.542800e+04)*expf((2.444000000000000e-01f*v_old_)))-(6.948000e-06*expf(((-4.391000000000000e-02f)*v_old_))))*(v_old_+3.778000000000000e+01f))/(1.000000000000000e+00f+expf((3.110000000000000e-01f*(v_old_+7.923000000000000e+01f)))));    }
#define IFNUMBER_5(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = ((6.000000000000000e-01f*expf((5.700000000000000e-02f*v_old_)))/(1.000000000000000e+00f+expf(((-1.000000000000000e-01f)*(v_old_+3.200000000000000e+01f)))));    }  else{ (name) = ((2.424000000000000e-02f*expf(((-1.052000000000000e-02f)*v_old_)))/(1.000000000000000e+00f+expf(((-1.378000000000000e-01f)*(v_old_+4.014000000000000e+01f)))));    }
#define IFNUMBER_6(name)if((celltype==1.000000000000000e+00f)) { (name) = (GNaL_b*6.000000000000000e-01f);    }  else{ (name) = GNaL_b;    }
#define IFNUMBER_7(name)if((celltype==1.000000000000000e+00f)) { (name) = (1.000000000000000e+00f-(9.500000000000000e-01f/(1.000000000000000e+00f+expf(((v_old_+EKshift+7.000000000000000e+01f)/5.000000000000000e+00f)))));    }  else{ (name) = 1.000000000000000e+00f;    }
#define IFNUMBER_8(name)if((celltype==1.000000000000000e+00f)) { (name) = (Gto_b*2.000000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Gto_b*2.000000000000000e+00f);    } else{ (name) = Gto_b;    }
#define IFNUMBER_9(name)if((v_old_>=3.149780000000000e+01f)) { (name) = 1.000000000000000e+00f;    }  else{ (name) = (1.076300000000000e+00f*expf(((-1.007000000000000e+00f)*expf(((-8.290000000000000e-02f)*v_old_)))));    }
#define IFNUMBER_10(name)if((celltype==1.000000000000000e+00f)) { (name) = (PCa_b*1.200000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (PCa_b*2.000000000000000e+00f);    } else{ (name) = PCa_b;    }
#define IFNUMBER_11(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKr_b*1.300000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (GKr_b*8.000000000000000e-01f);    } else{ (name) = GKr_b;    }
#define IFNUMBER_12(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKs_b*1.400000000000000e+00f);    }  else{ (name) = GKs_b;    }
#define IFNUMBER_13(name)if((celltype==1.000000000000000e+00f)) { (name) = (GK1_b*1.200000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (GK1_b*1.300000000000000e+00f);    } else{ (name) = GK1_b;    }
#define IFNUMBER_14(name)if((celltype==1.000000000000000e+00f)) { (name) = (Gncx_b*1.100000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Gncx_b*1.400000000000000e+00f);    } else{ (name) = Gncx_b;    }
#define IFNUMBER_15(name)if((celltype==1.000000000000000e+00f)) { (name) = (Pnak_b*9.000000000000000e-01f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Pnak_b*7.000000000000000e-01f);    } else{ (name) = Pnak_b;    }
#define IFNUMBER_16(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKb_b*6.000000000000000e-01f);    }  else{ (name) = GKb_b;    }
#define IFNUMBER_17(name)if((celltype==2.000000000000000e+00f)) { (name) = (calc_Jrel_inf_b*1.700000000000000e+00f);    }  else{ (name) = calc_Jrel_inf_b;    }
#define IFNUMBER_18(name)if((calc_tau_rel_b<1.000000000000000e-03f)) { (name) = 1.000000000000000e-03f;    }  else{ (name) = calc_tau_rel_b;    }
#define IFNUMBER_19(name)if((celltype==2.000000000000000e+00f)) { (name) = (calc_Jrel_infp_b*1.700000000000000e+00f);    }  else{ (name) = calc_Jrel_infp_b;    }
#define IFNUMBER_20(name)if((calc_tau_relp_b<1.000000000000000e-03f)) { (name) = 1.000000000000000e-03f;    }  else{ (name) = calc_tau_relp_b;    }
#define IFNUMBER_21(name)if((celltype==1.000000000000000e+00f)) { (name) = 1.300000000000000e+00f;    }  else{ (name) = 1.000000000000000e+00f;    }

#ifdef __CUDACC__

__constant__  size_t pitch;
__constant__  real abstol;
__constant__  real reltol;
__constant__  real max_dt;
__constant__  real min_dt;
__constant__  uint8_t use_adpt;

size_t pitch_h;

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes);

__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps);

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int thread_id, real dt);
inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt);
inline __device__ void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int thread_id);

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
#endif //MONOALG3D_MODEL_TORORD_FKATP_ENDO_H

