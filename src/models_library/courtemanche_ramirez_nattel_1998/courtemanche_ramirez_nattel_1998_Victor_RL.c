#include "courtemanche_ramirez_nattel_1998_Victor_RL.h"
#include <stdlib.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using courtemanche_ramirez_nattel_1998_Victor_RL CPU model\n");

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
        log_info("Using Rush-Larsen/Euler model to solve the ODEs\n");
    }

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {

        real *sv = &solver->sv[i * NEQ];
        if (solver->ode_extra_data == NULL) {
            sv[0]  = -8.118000e+01f; //V millivolt
            sv[1]  = 2.908000e-03f;  //m dimensionless
            sv[2]  = 9.649000e-01f;  //h dimensionless
            sv[3]  = 9.775000e-01f;  //j dimensionless
            sv[4]  = 3.043000e-02f;  //oa dimensionless
            sv[5]  = 9.992000e-01f;  //oi dimensionless
            sv[6]  = 4.966000e-03f;  //ua dimensionless
            sv[7]  = 9.986000e-01f;  //ui dimensionless
            sv[8]  = 3.296000e-05f;  //xr dimensionless
            sv[9]  = 1.869000e-02f;  //xs dimensionless
            sv[10] = 1.367000e-04f;  //d dimensionless
            sv[11] = 9.996000e-01f;  //f dimensionless
            sv[12] = 7.755000e-01f;  //f_Ca dimensionless
            sv[13] = 0.000000e+00f;  //u dimensionless
            sv[14] = 1.000000e+00f;  //v dimensionless
            sv[15] = 9.992000e-01f;  //w dimensionless
            sv[16] = 1.117000e+01f;  //Na_i millimolar
            sv[17] = 1.390000e+02f;  //K_i millimolar
            sv[18] = 1.013000e-04f;  //Ca_i millimolar
            sv[19] = 1.488000e+00f;  //Ca_up millimolar
            sv[20] = 1.488000e+00f;  //Ca_rel millimolar
        }
        else {
            real *extra_parameters = (real *)solver->ode_extra_data;
            sv[ 0] = extra_parameters[13];
            sv[ 1] = extra_parameters[14];
            sv[ 2] = extra_parameters[15];
            sv[ 3] = extra_parameters[16];
            sv[ 4] = extra_parameters[17];
            sv[ 5] = extra_parameters[18];
            sv[ 6] = extra_parameters[19];
            sv[ 7] = extra_parameters[20];
            sv[ 8] = extra_parameters[21];
            sv[ 9] = extra_parameters[22];
            sv[10] = extra_parameters[23];
            sv[11] = extra_parameters[24];
            sv[12] = extra_parameters[25];
            sv[13] = extra_parameters[26];
            sv[14] = extra_parameters[27];
            sv[15] = extra_parameters[28];
            sv[16] = extra_parameters[29];
            sv[17] = extra_parameters[30];
            sv[18] = extra_parameters[31];
            sv[19] = extra_parameters[32];
            sv[20] = extra_parameters[33];
        }
    }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    int offset = 34;
    bool adpt = ode_solver->adaptive;
    real *extra_parameters = NULL;

    if(ode_solver->ode_extra_data) {
        extra_parameters = ((real*)ode_solver->ode_extra_data);
    }
    else {
        log_error_and_exit("You need to specify a mask function when using this mixed model!\n");
    }

    #pragma omp parallel for private(sv_id)
    for (int i = 0; i < num_cells_to_solve; i++) {
        int mapping = extra_parameters[i+offset];

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver, ode_solver->ode_extra_data, mapping);
        }
        else {
            for (int j = 0; j < num_steps; ++j) {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], ode_solver->ode_extra_data, mapping);
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current, real *extra_parameters, int mapping) {

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, extra_parameters, mapping);

    // This model uses the Rush-Larsen method to solve the ODEs
    sv[0] = dt*rDY[0] + rY[0];          // Euler

    sv[1]  = rDY[1];                    // Rush-Larsen
    sv[2]  = rDY[2];                    // Rush-Larsen
    sv[3]  = rDY[3];                    // Rush-Larsen
    sv[4]  = rDY[4];                    // Rush-Larsen
    sv[5]  = rDY[5];                    // Rush-Larsen
    sv[6]  = rDY[6];                    // Rush-Larsen
    sv[7]  = rDY[7];                    // Rush-Larsen
    sv[8]  = rDY[8];                    // Rush-Larsen
    sv[9]  = rDY[9];                    // Rush-Larsen
    sv[10] = rDY[10];                   // Rush-Larsen
    sv[11] = rDY[11];                   // Rush-Larsen
    sv[12] = rDY[12];                   // Rush-Larsen
    sv[13] = rDY[13];                   // Rush-Larsen
    sv[14] = rDY[14];                   // Rush-Larsen
    sv[15] = rDY[15];                   // Rush-Larsen

    sv[16] = dt*rDY[16] + rY[16];       // Euler
    sv[17] = dt*rDY[17] + rY[17];       // Euler
    sv[18] = dt*rDY[18] + rY[18];       // Euler
    sv[19] = dt*rDY[19] + rY[19];       // Euler
    sv[20] = dt*rDY[20] + rY[20];       // Euler

}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver, real *extra_parameters, int mapping) {

    const real _beta_safety_ = 0.8;
    int numEDO = NEQ;

    real rDY[numEDO];

    real _tolerances_[numEDO];
    real _aux_tol = 0.0;
    // initializes the variables
    solver->ode_previous_dt[sv_id] = solver->ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k2__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k_aux__;

    real *dt = &solver->ode_dt[sv_id];
    real *time_new = &solver->ode_time_new[sv_id];
    real *previous_dt = &solver->ode_previous_dt[sv_id];

    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    RHS_cpu(sv, rDY, stim_curr, *dt, extra_parameters, mapping);
    *time_new += *dt;

    for(int i = 0; i < numEDO; i++) {
        _k1__[i] = rDY[i];
    }

    const real rel_tol = solver->rel_tol;
    const real abs_tol = solver->abs_tol;

    const real __tiny_ = pow(abs_tol, 2.0);

    real min_dt = solver->min_dt;
    real max_dt = solver->max_dt;

    while(1) {

        for(int i = 0; i < numEDO; i++) {
            // stores the old variables in a vector
            edos_old_aux_[i] = sv[i];
            // computes euler method
            edos_new_euler_[i] = _k1__[i] * *dt + edos_old_aux_[i];
            // steps ahead to compute the rk2 method
            sv[i] = edos_new_euler_[i];
        }

        *time_new += *dt;
        RHS_cpu(sv, rDY, stim_curr, *dt, extra_parameters, mapping);
        *time_new -= *dt; // step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < numEDO; i++) {
            // stores the new evaluation
            _k2__[i] = rDY[i];
            _aux_tol = fabs(edos_new_euler_[i]) * rel_tol;
            _tolerances_[i] = (abs_tol > _aux_tol) ? abs_tol : _aux_tol;
            // finds the greatest error between  the steps
            auxError = fabs(((*dt / 2.0) * (_k1__[i] - _k2__[i])) / _tolerances_[i]);

            greatestError = (auxError > greatestError) ? auxError : greatestError;
        }
        /// adapt the time step
        greatestError += __tiny_;
        *previous_dt = *dt;
        /// adapt the time step
        *dt = _beta_safety_ * (*dt) * sqrt(1.0f / greatestError);

        if(*dt < min_dt) {
            *dt = min_dt;
        } else if(*dt > max_dt) {
            *dt = max_dt;
        }

        if(*time_new + *dt > final_time) {
            *dt = final_time - *time_new;
        }

        // it doesn't accept the solution
        if(greatestError >= 1.0f && *dt > min_dt) {
            // restore the old values to do it again
            for(int i = 0; i < numEDO; i++) {
                sv[i] = edos_old_aux_[i];
            }
            // throw the results away and compute again
        } else {
            // it accepts the solutions
            if(greatestError >= 1.0) {
                printf("Accepting solution with error > %lf \n", greatestError);
            }

            _k_aux__ = _k2__;
            _k2__ = _k1__;
            _k1__ = _k_aux__;

            // it steps the method ahead, with euler solution
            for(int i = 0; i < numEDO; i++) {
                sv[i] = edos_new_euler_[i];
            }

            if(*time_new + *previous_dt >= final_time) {
                if(final_time == *time_new) {
                    break;
                } else if(*time_new < final_time) {
                    *dt = *previous_dt = final_time - *time_new;
                    *time_new += *previous_dt;
                    break;
                }
            } else {
                *time_new += *previous_dt;
            }
        }
    }

    free(_k1__);
    free(_k2__);
}

void RHS_cpu(const real *sv, real *rDY, real stim_current, real dt, real *extra_parameters, int mapping) {

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

    #include "courtemanche_ramirez_nattel_1998_Victor_RL_common.inc.c"
}


/*
[extra_data]
; Membrane model parameters
;GKur = 1.0
;GKr  = 1.0
;Gto  = 0.4
;GK1  = 2.0
;GKs  = 1.0
;GCaL = 0.35
;GNaK = 1.0
;GNCX = 1.0
;GNa  = 1.0
;GK2P = 0.0
;Gup  = 1.0
;Grel = 1.0
;Gleak= 1.0

; Initial conditions
;Vm = -81.2
;INa_va =    2.91e-03
;INa_vi_1 =  9.65e-01
;INa_vi_2 =  9.78e-01
;Ito_va =    3.04e-02
;Ito_vi =    9.99e-01
;ICaL_va =   1.37e-04
;ICaL_vi =   9.99e-01
;ICaL_ci =   7.75e-01
;IKur_va =   4.96e-03
;IKur_viS =  9.99e-01
;IKur_viF =  9.99e-01
;IKr_va =    3.29e-05
;IKs_va =    1.87e-02
;IK2P_va =   0.202895277583553
;CajSR =     1.49
;CanSR =     1.49
;RyRo =      1.06460418123529e-31
;RyRr =      1.0
;RyRi =      9.99e-01
;Nai =       1.12e+01
;Nao =       140
;Ki =        1.39e02
;Ko =        5.4
;Cai =       1.02e-04
;Cao =       1.8
*/