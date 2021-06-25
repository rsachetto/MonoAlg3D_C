#include "courtemanche_ramirez_nattel_1998.h"
#include <stdlib.h>

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using courtemanche_ramirez_nattel_1998 CPU model\n");

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

        sv[0] = -8.118000e+01f; //V millivolt
        sv[1] = 2.908000e-03f; //m dimensionless
        sv[2] = 9.649000e-01f; //h dimensionless
        sv[3] = 9.775000e-01f; //j dimensionless
        sv[4] = 3.043000e-02f; //oa dimensionless
        sv[5] = 9.992000e-01f; //oi dimensionless
        sv[6] = 4.966000e-03f; //ua dimensionless
        sv[7] = 9.986000e-01f; //ui dimensionless
        sv[8] = 3.296000e-05f; //xr dimensionless
        sv[9] = 1.869000e-02f; //xs dimensionless
        sv[10] = 1.367000e-04f; //d dimensionless
        sv[11] = 9.996000e-01f; //f dimensionless
        sv[12] = 7.755000e-01f; //f_Ca dimensionless
        sv[13] = 0.0; //u dimensionless
        sv[14] = 1.000000e+00f; //v dimensionless
        sv[15] = 9.992000e-01f; //w dimensionless
        sv[16] = 1.117000e+01f; //Na_i millimolar
        sv[17] = 1.390000e+02f; //K_i millimolar
        sv[18] = 1.013000e-04f; //Ca_i millimolar
        sv[19] = 1.488000e+00f; //Ca_up millimolar
        sv[20] = 1.488000e+00f; //Ca_rel millimolar
    }
}

void RHS_cpu(const real *sv, real *rDY, real stim_current, real dt) {

    //State variables
    const real V_old_ = sv[0];
    const real m_old_ = sv[1];
    const real h_old_ = sv[2];
    const real j_old_ = sv[3];
    const real oa_old_ = sv[4];
    const real oi_old_ = sv[5];
    const real ua_old_ = sv[6];
    const real ui_old_ = sv[7];
    const real xr_old_ = sv[8];
    const real xs_old_ = sv[9];
    const real d_old_ = sv[10];
    const real f_old_ = sv[11];
    const real f_Ca_old_ = sv[12];
    const real u_old_ = sv[13];
    const real v_old_ = sv[14];
    const real w_old_ = sv[15];
    const real Na_i_old_ = sv[16];
    const real K_i_old_ = sv[17];
    const real Ca_i_old_ = sv[18];
    const real Ca_up_old_ = sv[19];
    const real Ca_rel_old_ = sv[20];

    #include "courtemanche_ramirez_nattel_1998_common.inc.c"
}

#include "../default_solvers.c"
