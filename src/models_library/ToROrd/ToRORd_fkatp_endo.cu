#include "ToRORd_fkatp_endo.h"
#include <stddef.h>
#include <stdint.h>

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if(threadID < num_volumes) {

        *((real *)((char *)sv + pitch * 0) + threadID)  = -88.8691566357934;     // v millivolt
        *((real *)((char *)sv + pitch * 1) + threadID)  = 0.0110752904836162;    // CaMKt millimolar
        *((real *)((char *)sv + pitch * 2) + threadID)  = 12.0996647655188;      // nai millimolar
        *((real *)((char *)sv + pitch * 3) + threadID)  = 12.1000028563765;      // nass millimolar
        *((real *)((char *)sv + pitch * 4) + threadID)  = 142.412524737626;      // ki millimolar
        *((real *)((char *)sv + pitch * 5) + threadID)  = 142.412481425842;      // kss millimolar
        *((real *)((char *)sv + pitch * 6) + threadID)  = 7.45541572746214e-05;  // cai millimolar
        *((real *)((char *)sv + pitch * 7) + threadID)  = 6.50418928341426e-05;  // cass millimolar
        *((real *)((char *)sv + pitch * 8) + threadID)  = 1.53037019085812;      // cansr millimolar
        *((real *)((char *)sv + pitch * 9) + threadID)  = 1.52803094224238;      // cajsr millimolar
        *((real *)((char *)sv + pitch * 10) + threadID) = 0.000787657400526199;  // m dimensionless
        *((real *)((char *)sv + pitch * 11) + threadID) = 0.830658198588696;     // h dimensionless
        *((real *)((char *)sv + pitch * 12) + threadID) = 0.830466744399495;     // j dimensionless
        *((real *)((char *)sv + pitch * 13) + threadID) = 0.674096901201792;     // hp dimensionless
        *((real *)((char *)sv + pitch * 14) + threadID) = 0.830093612199637;     // jp dimensionless
        *((real *)((char *)sv + pitch * 15) + threadID) = 0.000159670117055769;  // mL dimensionless
        *((real *)((char *)sv + pitch * 16) + threadID) = 0.528261721740178;     // hL dimensionless
        *((real *)((char *)sv + pitch * 17) + threadID) = 0.288775833197764;     // hLp dimensionless
        *((real *)((char *)sv + pitch * 18) + threadID) = 0.000944249645410894;  // a dimensionless
        *((real *)((char *)sv + pitch * 19) + threadID) = 0.999616956857814;     // iF dimensionless
        *((real *)((char *)sv + pitch * 20) + threadID) = 0.593680589620082;     // iS dimensionless
        *((real *)((char *)sv + pitch * 21) + threadID) = 0.000481107253796778;  // ap dimensionless
        *((real *)((char *)sv + pitch * 22) + threadID) = 0.999616964658062;     // iFp dimensionless
        *((real *)((char *)sv + pitch * 23) + threadID) = 0.654092074678260;     // iSp dimensionless
        *((real *)((char *)sv + pitch * 24) + threadID) = 8.86091322819384e-29;  // d dimensionless
        *((real *)((char *)sv + pitch * 25) + threadID) = 0.999999992783113;     // ff dimensionless
        *((real *)((char *)sv + pitch * 26) + threadID) = 0.938965241412012;     // fs dimensionless
        *((real *)((char *)sv + pitch * 27) + threadID) = 0.999999992783179;     // fcaf dimensionless
        *((real *)((char *)sv + pitch * 28) + threadID) = 0.999900458262832;     // fcas dimensionless
        *((real *)((char *)sv + pitch * 29) + threadID) = 0.999977476316330;     // jca dimensionless
        *((real *)((char *)sv + pitch * 30) + threadID) = 0.999999992566681;     // ffp dimensionless
        *((real *)((char *)sv + pitch * 31) + threadID) = 0.999999992766279;     // fcafp dimensionless
        *((real *)((char *)sv + pitch * 32) + threadID) = 0.000492094765239740;  // nca_ss dimensionless
        *((real *)((char *)sv + pitch * 33) + threadID) = 0.000833711885764158;  // nca_i dimensionless
        *((real *)((char *)sv + pitch * 34) + threadID) = 0.998073652444028;     // C3 dimensionless
        *((real *)((char *)sv + pitch * 35) + threadID) = 0.000844745297078649;  // C2 dimensionless
        *((real *)((char *)sv + pitch * 36) + threadID) = 0.000698171876592920;  // C1 dimensionless
        *((real *)((char *)sv + pitch * 37) + threadID) = 0.000370404872169913;  // O dimensionless
        *((real *)((char *)sv + pitch * 38) + threadID) = 1.30239063420973e-05;  // I dimensionless
        *((real *)((char *)sv + pitch * 39) + threadID) = 0.247156543918935;     // xs1 dimensionless
        *((real *)((char *)sv + pitch * 40) + threadID) = 0.000175017075236424;  // xs2 dimensionless
        *((real *)((char *)sv + pitch * 41) + threadID) = 3.90843796133124e-24;  // Jrel_np millimolar_per_millisecond
        *((real *)((char *)sv + pitch * 42) + threadID) = -1.88428892080206e-22; // Jrel_p millimolar_per_millisecond

        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * NEQ) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * (NEQ + 1)) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * (NEQ + 2)) + threadID) = 0.0;    // previous dt
        }
    }
}

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt, size_t pitch, bool use_adpt) {

    // State variables
    real v_old_;
    real CaMKt_old_;
    real nai_old_;
    real nass_old_;
    real ki_old_;
    real kss_old_;
    real cai_old_;
    real cass_old_;
    real cansr_old_;
    real cajsr_old_;
    real m_old_;
    real h_old_;
    real j_old_;
    real hp_old_;
    real jp_old_;
    real mL_old_;
    real hL_old_;
    real hLp_old_;
    real a_old_;
    real iF_old_;
    real iS_old_;
    real ap_old_;
    real iFp_old_;
    real iSp_old_;
    real d_old_;
    real ff_old_;
    real fs_old_;
    real fcaf_old_;
    real fcas_old_;
    real jca_old_;
    real ffp_old_;
    real fcafp_old_;
    real nca_ss_old_;
    real nca_i_old_;
    real C3_old_;
    real C2_old_;
    real C1_old_;
    real O_old_;
    real I_old_;
    real xs1_old_;
    real xs2_old_;
    real Jrel_np_old_;
    real Jrel_p_old_;

    if(use_adpt) {
        v_old_ = sv[0];
        CaMKt_old_ = sv[1];
        nai_old_ = sv[2];
        nass_old_ = sv[3];
        ki_old_ = sv[4];
        kss_old_ = sv[5];
        cai_old_ = sv[6];
        cass_old_ = sv[7];
        cansr_old_ = sv[8];
        cajsr_old_ = sv[9];
        m_old_ = sv[10];
        h_old_ = sv[11];
        j_old_ = sv[12];
        hp_old_ = sv[13];
        jp_old_ = sv[14];
        mL_old_ = sv[15];
        hL_old_ = sv[16];
        hLp_old_ = sv[17];
        a_old_ = sv[18];
        iF_old_ = sv[19];
        iS_old_ = sv[20];
        ap_old_ = sv[21];
        iFp_old_ = sv[22];
        iSp_old_ = sv[23];
        d_old_ = sv[24];
        ff_old_ = sv[25];
        fs_old_ = sv[26];
        fcaf_old_ = sv[27];
        fcas_old_ = sv[28];
        jca_old_ = sv[29];
        ffp_old_ = sv[30];
        fcafp_old_ = sv[31];
        nca_ss_old_ = sv[32];
        nca_i_old_ = sv[33];
        C3_old_ = sv[34];
        C2_old_ = sv[35];
        C1_old_ = sv[36];
        O_old_ = sv[37];
        I_old_ = sv[38];
        xs1_old_ = sv[39];
        xs2_old_ = sv[40];
        Jrel_np_old_ = sv[41];
        Jrel_p_old_ = sv[42];
    } else {
        //    //State variables
        v_old_ = *((real *)((char *)sv + pitch * 0) + threadID_);
        CaMKt_old_ = *((real *)((char *)sv + pitch * 1) + threadID_);
        nai_old_ = *((real *)((char *)sv + pitch * 2) + threadID_);
        nass_old_ = *((real *)((char *)sv + pitch * 3) + threadID_);
        ki_old_ = *((real *)((char *)sv + pitch * 4) + threadID_);
        kss_old_ = *((real *)((char *)sv + pitch * 5) + threadID_);
        cai_old_ = *((real *)((char *)sv + pitch * 6) + threadID_);
        cass_old_ = *((real *)((char *)sv + pitch * 7) + threadID_);
        cansr_old_ = *((real *)((char *)sv + pitch * 8) + threadID_);
        cajsr_old_ = *((real *)((char *)sv + pitch * 9) + threadID_);
        m_old_ = *((real *)((char *)sv + pitch * 10) + threadID_);
        h_old_ = *((real *)((char *)sv + pitch * 11) + threadID_);
        j_old_ = *((real *)((char *)sv + pitch * 12) + threadID_);
        hp_old_ = *((real *)((char *)sv + pitch * 13) + threadID_);
        jp_old_ = *((real *)((char *)sv + pitch * 14) + threadID_);
        mL_old_ = *((real *)((char *)sv + pitch * 15) + threadID_);
        hL_old_ = *((real *)((char *)sv + pitch * 16) + threadID_);
        hLp_old_ = *((real *)((char *)sv + pitch * 17) + threadID_);
        a_old_ = *((real *)((char *)sv + pitch * 18) + threadID_);
        iF_old_ = *((real *)((char *)sv + pitch * 19) + threadID_);
        iS_old_ = *((real *)((char *)sv + pitch * 20) + threadID_);
        ap_old_ = *((real *)((char *)sv + pitch * 21) + threadID_);
        iFp_old_ = *((real *)((char *)sv + pitch * 22) + threadID_);
        iSp_old_ = *((real *)((char *)sv + pitch * 23) + threadID_);
        d_old_ = *((real *)((char *)sv + pitch * 24) + threadID_);
        ff_old_ = *((real *)((char *)sv + pitch * 25) + threadID_);
        fs_old_ = *((real *)((char *)sv + pitch * 26) + threadID_);
        fcaf_old_ = *((real *)((char *)sv + pitch * 27) + threadID_);
        fcas_old_ = *((real *)((char *)sv + pitch * 28) + threadID_);
        jca_old_ = *((real *)((char *)sv + pitch * 29) + threadID_);
        ffp_old_ = *((real *)((char *)sv + pitch * 30) + threadID_);
        fcafp_old_ = *((real *)((char *)sv + pitch * 31) + threadID_);
        nca_ss_old_ = *((real *)((char *)sv + pitch * 32) + threadID_);
        nca_i_old_ = *((real *)((char *)sv + pitch * 33) + threadID_);
        C3_old_ = *((real *)((char *)sv + pitch * 34) + threadID_);
        C2_old_ = *((real *)((char *)sv + pitch * 35) + threadID_);
        C1_old_ = *((real *)((char *)sv + pitch * 36) + threadID_);
        O_old_ = *((real *)((char *)sv + pitch * 37) + threadID_);
        I_old_ = *((real *)((char *)sv + pitch * 38) + threadID_);
        xs1_old_ = *((real *)((char *)sv + pitch * 39) + threadID_);
        xs2_old_ = *((real *)((char *)sv + pitch * 40) + threadID_);
        Jrel_np_old_ = *((real *)((char *)sv + pitch * 41) + threadID_);
        Jrel_p_old_ = *((real *)((char *)sv + pitch * 42) + threadID_);
    }

    #include "ToROrd_common.inc.c"
}

//Include the default solver used by all models.
#include "../default_solvers.cu"
