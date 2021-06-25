#include "ToRORd_fkatp_endo.h"
#include <stdlib.h>

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using ToRORd_fkatp_endo CPU model\n");

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    bool adpt = solver->adaptive;

    if(adpt) {
        solver->ode_dt = (real*)malloc(num_cells*sizeof(real));

        OMP(parallel for)
        for(int i = 0; i < num_cells; i++) {
            solver->ode_dt[i] = solver->min_dt;
        }

        solver->ode_previous_dt = (real*)calloc(num_cells, sizeof(real));
        solver->ode_time_new    = (real*)calloc(num_cells, sizeof(real));
        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_info("Using Euler model to solve the ODEs\n");
    }


    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {

        real *sv = &solver->sv[i * NEQ];

        sv[0] = -8.876380e+01f; //v millivolt
        sv[1] = 1.110000e-02f; //CaMKt millimolar
        sv[2] = 1.210250e+01f; //nai millimolar
        sv[3] = 1.210290e+01f; //nass millimolar
        sv[4] = 1.423002e+02f; //ki millimolar
        sv[5] = 1.423002e+02f; //kss millimolar
        sv[6] = 8.158300e-05f; //cai millimolar
        sv[7] = 7.030500e-05f; //cass millimolar
        sv[8] = 1.521100e+00f; //cansr millimolar
        sv[9] = 1.521400e+00f; //cajsr millimolar
        sv[10] = 8.057200e-04f; //m dimensionless
        sv[11] = 8.286000e-01f; //h dimensionless
        sv[12] = 8.284000e-01f; //j dimensionless
        sv[13] = 6.707000e-01f; //hp dimensionless
        sv[14] = 8.281000e-01f; //jp dimensionless
        sv[15] = 1.629000e-04f; //mL dimensionless
        sv[16] = 5.255000e-01f; //hL dimensionless
        sv[17] = 2.872000e-01f; //hLp dimensionless
        sv[18] = 9.509800e-04f; //a dimensionless
        sv[19] = 9.996000e-01f; //iF dimensionless
        sv[20] = 5.936000e-01f; //iS dimensionless
        sv[21] = 4.845400e-04f; //ap dimensionless
        sv[22] = 9.996000e-01f; //iFp dimensionless
        sv[23] = 6.538000e-01f; //iSp dimensionless
        sv[24] = 8.108400e-09f; //d dimensionless
        sv[25] = 1.000000e+00f; //ff dimensionless
        sv[26] = 9.390000e-01f; //fs dimensionless
        sv[27] = 1.000000e+00f; //fcaf dimensionless
        sv[28] = 9.999000e-01f; //fcas dimensionless
        sv[29] = 1.000000e+00f; //jca dimensionless
        sv[30] = 1.000000e+00f; //ffp dimensionless
        sv[31] = 1.000000e+00f; //fcafp dimensionless
        sv[32] = 6.646200e-04f; //nca_ss dimensionless
        sv[33] = 1.200000e-03f; //nca_i dimensionless
        sv[34] = 9.981000e-01f; //C3 dimensionless
        sv[35] = 8.510900e-04f; //C2 dimensionless
        sv[36] = 7.034400e-04f; //C1 dimensionless
        sv[37] = 3.758500e-04f; //O dimensionless
        sv[38] = 1.328900e-05f; //I dimensionless
        sv[39] = 2.480000e-01f; //xs1 dimensionless
        sv[40] = 1.770700e-04f; //xs2 dimensionless
        sv[41] = 1.612900e-22f; //Jrel_np millimolar_per_millisecond
        sv[42] = 1.247500e-20f; //Jrel_p millimolar_per_millisecond
    }
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {

    //State variables
    const real v_old_ = sv[0];
    const real CaMKt_old_ = sv[1];
    const real nai_old_ = sv[2];
    const real nass_old_ = sv[3];
    const real ki_old_ = sv[4];
    const real kss_old_ = sv[5];
    const real cai_old_ = sv[6];
    const real cass_old_ = sv[7];
    const real cansr_old_ = sv[8];
    const real cajsr_old_ = sv[9];
    const real m_old_ = sv[10];
    const real h_old_ = sv[11];
    const real j_old_ = sv[12];
    const real hp_old_ = sv[13];
    const real jp_old_ = sv[14];
    const real mL_old_ = sv[15];
    const real hL_old_ = sv[16];
    const real hLp_old_ = sv[17];
    const real a_old_ = sv[18];
    const real iF_old_ = sv[19];
    const real iS_old_ = sv[20];
    const real ap_old_ = sv[21];
    const real iFp_old_ = sv[22];
    const real iSp_old_ = sv[23];
    const real d_old_ = sv[24];
    const real ff_old_ = sv[25];
    const real fs_old_ = sv[26];
    const real fcaf_old_ = sv[27];
    const real fcas_old_ = sv[28];
    const real jca_old_ = sv[29];
    const real ffp_old_ = sv[30];
    const real fcafp_old_ = sv[31];
    const real nca_ss_old_ = sv[32];
    const real nca_i_old_ = sv[33];
    const real C3_old_ = sv[34];
    const real C2_old_ = sv[35];
    const real C1_old_ = sv[36];
    const real O_old_ = sv[37];
    const real I_old_ = sv[38];
    const real xs1_old_ = sv[39];
    const real xs2_old_ = sv[40];
    const real Jrel_np_old_ = sv[41];
    const real Jrel_p_old_ = sv[42];

    #include "ToROrd_common.inc.c"
}

#include "../default_solvers.c"
