#include "courtemanche_ramirez_nattel_1998.h"
#include <stddef.h>
#include <stdint.h>

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

         *((real * )((char *) sv + pitch * 0) + threadID) = -8.118000e+01f; //V millivolt 
         *((real * )((char *) sv + pitch * 1) + threadID) = 2.908000e-03f; //m dimensionless 
         *((real * )((char *) sv + pitch * 2) + threadID) = 9.649000e-01f; //h dimensionless 
         *((real * )((char *) sv + pitch * 3) + threadID) = 9.775000e-01f; //j dimensionless 
         *((real * )((char *) sv + pitch * 4) + threadID) = 3.043000e-02f; //oa dimensionless 
         *((real * )((char *) sv + pitch * 5) + threadID) = 9.992000e-01f; //oi dimensionless 
         *((real * )((char *) sv + pitch * 6) + threadID) = 4.966000e-03f; //ua dimensionless 
         *((real * )((char *) sv + pitch * 7) + threadID) = 9.986000e-01f; //ui dimensionless 
         *((real * )((char *) sv + pitch * 8) + threadID) = 3.296000e-05f; //xr dimensionless 
         *((real * )((char *) sv + pitch * 9) + threadID) = 1.869000e-02f; //xs dimensionless 
         *((real * )((char *) sv + pitch * 10) + threadID) = 1.367000e-04f; //d dimensionless 
         *((real * )((char *) sv + pitch * 11) + threadID) = 9.996000e-01f; //f dimensionless 
         *((real * )((char *) sv + pitch * 12) + threadID) = 7.755000e-01f; //f_Ca dimensionless 
         *((real * )((char *) sv + pitch * 13) + threadID) = 0.0f; //u dimensionless
         *((real * )((char *) sv + pitch * 14) + threadID) = 1.000000e+00f; //v dimensionless 
         *((real * )((char *) sv + pitch * 15) + threadID) = 9.992000e-01f; //w dimensionless 
         *((real * )((char *) sv + pitch * 16) + threadID) = 1.117000e+01f; //Na_i millimolar 
         *((real * )((char *) sv + pitch * 17) + threadID) = 1.390000e+02f; //K_i millimolar 
         *((real * )((char *) sv + pitch * 18) + threadID) = 1.013000e-04f; //Ca_i millimolar 
         *((real * )((char *) sv + pitch * 19) + threadID) = 1.488000e+00f; //Ca_up millimolar 
         *((real * )((char *) sv + pitch * 20) + threadID) = 1.488000e+00f; //Ca_rel millimolar 

        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * NEQ) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * (NEQ + 1)) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * (NEQ + 2)) + threadID) = 0.0;    // previous dt
        }
	}
}

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int thread_id, real dt, size_t pitch, bool use_adpt_dt) {

    //State variables
    real V_old_; //millivolt
    real m_old_; //dimensionless
    real h_old_; //dimensionless
    real j_old_; //dimensionless
    real oa_old_; //dimensionless
    real oi_old_; //dimensionless
    real ua_old_; //dimensionless
    real ui_old_; //dimensionless
    real xr_old_; //dimensionless
    real xs_old_; //dimensionless
    real d_old_; //dimensionless
    real f_old_; //dimensionless
    real f_Ca_old_; //dimensionless
    real u_old_; //dimensionless
    real v_old_; //dimensionless
    real w_old_; //dimensionless
    real Na_i_old_; //millimolar
    real K_i_old_; //millimolar
    real Ca_i_old_; //millimolar
    real Ca_up_old_; //millimolar
    real Ca_rel_old_; //millimolar


    if(use_adpt_dt) {
        V_old_ =  sv[0];
        m_old_ =  sv[1];
        h_old_ =  sv[2];
        j_old_ =  sv[3];
        oa_old_ =  sv[4];
        oi_old_ =  sv[5];
        ua_old_ =  sv[6];
        ui_old_ =  sv[7];
        xr_old_ =  sv[8];
        xs_old_ =  sv[9];
        d_old_ =  sv[10];
        f_old_ =  sv[11];
        f_Ca_old_ =  sv[12];
        u_old_ =  sv[13];
        v_old_ =  sv[14];
        w_old_ =  sv[15];
        Na_i_old_ =  sv[16];
        K_i_old_ =  sv[17];
        Ca_i_old_ =  sv[18];
        Ca_up_old_ =  sv[19];
        Ca_rel_old_ =  sv[20];
    } else {
        V_old_ =  *((real*)((char*)sv + pitch * 0) + thread_id);
        m_old_ =  *((real*)((char*)sv + pitch * 1) + thread_id);
        h_old_ =  *((real*)((char*)sv + pitch * 2) + thread_id);
        j_old_ =  *((real*)((char*)sv + pitch * 3) + thread_id);
        oa_old_ =  *((real*)((char*)sv + pitch * 4) + thread_id);
        oi_old_ =  *((real*)((char*)sv + pitch * 5) + thread_id);
        ua_old_ =  *((real*)((char*)sv + pitch * 6) + thread_id);
        ui_old_ =  *((real*)((char*)sv + pitch * 7) + thread_id);
        xr_old_ =  *((real*)((char*)sv + pitch * 8) + thread_id);
        xs_old_ =  *((real*)((char*)sv + pitch * 9) + thread_id);
        d_old_ =  *((real*)((char*)sv + pitch * 10) + thread_id);
        f_old_ =  *((real*)((char*)sv + pitch * 11) + thread_id);
        f_Ca_old_ =  *((real*)((char*)sv + pitch * 12) + thread_id);
        u_old_ =  *((real*)((char*)sv + pitch * 13) + thread_id);
        v_old_ =  *((real*)((char*)sv + pitch * 14) + thread_id);
        w_old_ =  *((real*)((char*)sv + pitch * 15) + thread_id);
        Na_i_old_ =  *((real*)((char*)sv + pitch * 16) + thread_id);
        K_i_old_ =  *((real*)((char*)sv + pitch * 17) + thread_id);
        Ca_i_old_ =  *((real*)((char*)sv + pitch * 18) + thread_id);
        Ca_up_old_ =  *((real*)((char*)sv + pitch * 19) + thread_id);
        Ca_rel_old_ =  *((real*)((char*)sv + pitch * 20) + thread_id);
    }

    #include "courtemanche_ramirez_nattel_1998_common.inc.c"

}

//Include the default solver used by all models.
#include "../default_solvers.cu"
