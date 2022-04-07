#include "ToRORd_fkatp_mixed_endo_mid_epi.h"
#include <stdlib.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using ToRORd_fkatp_2019 CPU model\n");

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
        log_info("Using Adaptive timestep model to solve the ODEs\n");
    } else {
        log_info("Using Fixed timestep to solve the ODEs\n");
    }

    real *initial_endo = NULL;
    real *initial_epi = NULL;
    real *initial_mid = NULL;
    real *mapping = NULL;
    if(solver->ode_extra_data) {
        initial_endo = (real *)solver->ode_extra_data;
        initial_epi = (real *)solver->ode_extra_data+NEQ;
        initial_mid = (real *)solver->ode_extra_data+NEQ+NEQ;
        mapping = (real *)solver->ode_extra_data+NEQ+NEQ+NEQ;

        OMP(parallel for)
        for(uint32_t i = 0; i < num_cells; i++){
            
            real *sv = &solver->sv[i * NEQ];

            for (int j = 0; j < NEQ; j++) {
                if (mapping[i] == ENDO)
                    sv[j] = initial_endo[j];
                else if (mapping[i] == EPI)
                    sv[j] = initial_epi[j];
                else
                    sv[j] = initial_mid[j];
            }
        }
    }
    else {
        log_info("[INFO] You should supply a mask function to tag the cells when using this mixed model!\n");
        log_info("[INFO] Considering all cells ENDO!\n");
        
        OMP(parallel for)
        for(uint32_t i = 0; i < num_cells; i++){
            
            real *sv = &solver->sv[i * NEQ];

            // Default initial conditions (endocardium cell)
            /*
            sv[0] = -88.7638;
            sv[1] = 0.0111;
            sv[2] = 7.0305e-5;
            sv[3] = 12.1025;
            sv[4] = 12.1029;
            sv[5] = 142.3002;
            sv[6] = 142.3002;
            sv[7] = 1.5211;
            sv[8] = 1.5214;
            sv[9] = 8.1583e-05;
            sv[10] = 8.0572e-4;
            sv[11] = 0.8286;
            sv[12] = 0.8284;
            sv[13] = 0.6707;
            sv[14] = 0.8281;
            sv[15] = 1.629e-4;
            sv[16] = 0.5255;
            sv[17] = 0.2872;
            sv[18] = 9.5098e-4;
            sv[19] = 0.9996;
            sv[20] = 0.5936;
            sv[21] = 4.8454e-4;
            sv[22] = 0.9996;
            sv[23] = 0.6538;
            sv[24] = 8.1084e-9;
            sv[25] = 1.0;
            sv[26] = 0.939;
            sv[27] = 1.0;
            sv[28] = 0.9999;
            sv[29] = 1.0;
            sv[30] = 1.0;
            sv[31] = 1.0;
            sv[32] = 6.6462e-4;
            sv[33] = 0.0012;
            sv[34] = 7.0344e-4;
            sv[35] = 8.5109e-4;
            sv[36] = 0.9981;
            sv[37] = 1.3289e-5;
            sv[38] = 3.7585e-4;
            sv[39] = 0.248;
            sv[40] = 1.7707e-4;
            sv[41] = 1.6129e-22;
            sv[42] = 1.2475e-20;
            */

            // Steady-state after 200 beats (endocardium cell)
            sv[0] = -8.890585e+01;
            sv[1] = 1.107642e-02;
            sv[2] = 6.504164e-05;
            sv[3] = 1.210818e+01;
            sv[4] = 1.210851e+01;
            sv[5] = 1.426206e+02;
            sv[6] = 1.426205e+02;
            sv[7] = 1.530373e+00;
            sv[8] = 1.528032e+00;
            sv[9] = 7.455488e-05;
            sv[10] = 7.814592e-04;
            sv[11] = 8.313839e-01;
            sv[12] = 8.311938e-01;
            sv[13] = 6.752873e-01;
            sv[14] = 8.308255e-01;
            sv[15] = 1.585610e-04;
            sv[16] = 5.294475e-01;
            sv[17] = 2.896996e-01;
            sv[18] = 9.419166e-04;
            sv[19] = 9.996194e-01;
            sv[20] = 5.938602e-01;
            sv[21] = 4.799180e-04;
            sv[22] = 9.996194e-01;
            sv[23] = 6.543754e-01;
            sv[24] = -2.898677e-33;
            sv[25] = 1.000000e+00;
            sv[26] = 9.389659e-01;
            sv[27] = 1.000000e+00;
            sv[28] = 9.999003e-01;
            sv[29] = 9.999773e-01;
            sv[30] = 1.000000e+00;
            sv[31] = 1.000000e+00;
            sv[32] = 4.920606e-04;
            sv[33] = 8.337021e-04;
            sv[34] = 6.962775e-04;
            sv[35] = 8.425453e-04;
            sv[36] = 9.980807e-01;
            sv[37] = 1.289824e-05;
            sv[38] = 3.675442e-04;
            sv[39] = 2.471690e-01;
            sv[40] = 1.742987e-04;
            sv[41] = 5.421027e-24;
            sv[42] = 6.407933e-23;
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
    bool adpt = ode_solver->adaptive;

    // Get the mapping array
    real *mapping = NULL;
    if (ode_solver->ode_extra_data) {
        mapping = (real *)ode_solver->ode_extra_data+NEQ+NEQ+NEQ;
    }

    OMP(parallel for private(sv_id))
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            if (ode_solver->ode_extra_data) {
                //solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], mapping[i], current_t + dt, sv_id, ode_solver);
                solve_rush_larsen_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], mapping[i], current_t + dt, sv_id, ode_solver);
            }
            else {
                //solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], 0.0, current_t + dt, sv_id, ode_solver);
                solve_rush_larsen_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], 0.0, current_t + dt, sv_id, ode_solver);
            }
        }
        else {
            for (int j = 0; j < num_steps; ++j) {
                if (ode_solver->ode_extra_data) {
                    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], mapping[i]);
                }
                else {
                    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], 0.0);
                }
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current, real mapping) {

    // -------------------------------------------------------------------------------------------
    // MODEL SPECIFIC:
    // set the variables which are non-linear and hodkin-huxley type
    const real TOLERANCE = 1e-8;
    bool is_rush_larsen[NEQ];
    for (int i = 0; i < NEQ; i++) {
        is_rush_larsen[i] = ((i >= 10 && i <= 31) || (i >= 39 && i <= 42)) ? true : false;        
    }
    // -------------------------------------------------------------------------------------------
    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    // Compute 'a', 'b' coefficients alongside 'rhs'
    real a[NEQ], b[NEQ];
    RHS_RL_cpu(a, b, sv, rDY, stim_current, dt, mapping);

    // Solve variables based on its type:
    //  Non-linear = Euler
    //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
    for (int i = 0; i < NEQ; i++) {
        if (is_rush_larsen[i]) {
            if (abs(a[i]) < TOLERANCE) { 
                sv[i] = rY[i] + dt * (rY[i] * a[i] + b[i]);
            } 
            else {
                real aux = b[i] / a[i];
                sv[i] = exp(a[i] * dt)*(rY[i] + aux) - aux;
            }
        }
        else {
            sv[i] = dt * rDY[i] + rY[i];
        }
    }
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int sv_id, struct ode_solver *solver) {

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

    // Keep 'dt' inside the adaptive interval
    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    RHS_cpu(sv, rDY, stim_curr, *dt, mapping);
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
        RHS_cpu(sv, rDY, stim_curr, *dt, mapping);
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

void solve_rush_larsen_cpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int sv_id, struct ode_solver *solver) {
    
    // TODO: Remove this boolean array and write the full algebraics instead ...
    // -------------------------------------------------------------------------------------------
    // MODEL SPECIFIC:
    // set the variables which are non-linear and hodkin-huxley type
    const real TOLERANCE = 1e-09;
    bool is_rush_larsen[NEQ];
    for (int i = 0; i < NEQ; i++) {
        is_rush_larsen[i] = ((i >= 10 && i <= 31) || (i >= 39 && i <= 42)) ? true : false;        
    }
    // -------------------------------------------------------------------------------------------

    const real _beta_safety_ = 0.8;
    int numEDO = NEQ;

    real rDY[numEDO];

    // initializes the variables
    solver->ode_previous_dt[sv_id] = solver->ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k2__ = (real *)malloc(sizeof(real) * numEDO);
    real *a_ = (real *)malloc(sizeof(real) * numEDO);
    real *b_ = (real *)malloc(sizeof(real) * numEDO);
    real *a_new = (real *)malloc(sizeof(real) * numEDO);
    real *b_new = (real *)malloc(sizeof(real) * numEDO);
    real *_k_aux__, *_a_aux__, *_b_aux__;

    real *dt = &solver->ode_dt[sv_id];
    real *time_new = &solver->ode_time_new[sv_id];
    real *previous_dt = &solver->ode_previous_dt[sv_id];

    // Keep 'dt' inside the adaptive interval
    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    RHS_RL_cpu(a_, b_, sv, rDY, stim_curr, *dt, mapping);
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
            // computes euler/rush-larsen method
            if (is_rush_larsen[i])
                edos_new_euler_[i] = (fabs(a_[i]) < abs_tol) ? edos_old_aux_[i] + (edos_old_aux_[i] * a_[i] + b_[i])*(*dt) : \
                                                  exp(a_[i]*(*dt))*(edos_old_aux_[i] + (b_[i] / a_[i])) - (b_[i] / a_[i]);
            else
                edos_new_euler_[i] = _k1__[i] * *dt + edos_old_aux_[i];
            // steps ahead to compute the rk2 method
            sv[i] = edos_new_euler_[i];
        }

        *time_new += *dt;
        RHS_RL_cpu(a_new, b_new, sv, rDY, stim_curr, *dt, mapping);
        *time_new -= *dt; // step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < numEDO; i++) {
            _k2__[i] = rDY[i];
            if (is_rush_larsen[i]) {
                real as = (a_[i] + a_new[i]) * 0.5;
                real bs = (b_[i] + b_new[i]) * 0.5;
                real y_2nd_order = (fabs(as) < abs_tol) ? edos_old_aux_[i] + (*dt) * (edos_old_aux_[i]*as + bs) : \
                                                       exp(as*(*dt))*(edos_old_aux_[i] + (bs/as)) - (bs/as);
                auxError = (fabs(y_2nd_order) < abs_tol) ? fabs(edos_new_euler_[i] - abs_tol) : \
                                                        fabs( (y_2nd_order - edos_new_euler_[i])/(y_2nd_order) );
                greatestError = (auxError > greatestError) ? auxError : greatestError;
            }
            else {
                real f = (_k1__[i] + _k2__[i]) * 0.5;
                real y_2nd_order = edos_old_aux_[i] + (*dt) * f;
                auxError = (fabs(y_2nd_order) < abs_tol) ? fabs(edos_new_euler_[i] - abs_tol) : \
                                                        fabs( (y_2nd_order - edos_new_euler_[i])/(y_2nd_order) );
                greatestError = (auxError > greatestError) ? auxError : greatestError;
            }
        }
        /// adapt the time step
        greatestError += __tiny_;
        *previous_dt = *dt;
        /// adapt the time step
        *dt = (*dt) * sqrt(0.5 * rel_tol / greatestError);            // Jhonny`s formula

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

            // Swap pointers
            _k_aux__ = _k2__;
            _k2__ = _k1__;
            _k1__ = _k_aux__;

            _a_aux__ = a_;
            a_ = a_new;
            a_new = _a_aux__;

            _b_aux__ = b_;
            b_ = b_new;
            b_new = _b_aux__;

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
    free(a_);
    free(b_);
    free(a_new);
    free(b_new);
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real mapping) {

    // Get the celltype for the current cell
    real celltype = mapping;

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real v = sv[0];
    real CaMKt = sv[1];
    real cass = sv[2];
    real nai = sv[3];
    real nass = sv[4];
    real ki = sv[5];
    real kss = sv[6];
    real cansr = sv[7];
    real cajsr = sv[8];
    real cai = sv[9];
    real m = sv[10];
    real h = sv[11];
    real j = sv[12];
    real hp = sv[13];
    real jp = sv[14];
    real mL = sv[15];
    real hL = sv[16];
    real hLp = sv[17];
    real a = sv[18];
    real iF = sv[19];
    real iS = sv[20];
    real ap = sv[21];
    real iFp = sv[22];
    real iSp = sv[23];
    real d = sv[24];
    real ff = sv[25];
    real fs = sv[26];
    real fcaf = sv[27];
    real fcas = sv[28];
    real jca = sv[29];
    real ffp = sv[30];
    real fcafp = sv[31];
    real nca = sv[32];
    real nca_i = sv[33];
    real ikr_c0 = sv[34];
    real ikr_c1 = sv[35];
    real ikr_c2 = sv[36];
    real ikr_i = sv[37];
    real ikr_o = sv[38];
    real xs1 = sv[39];
    real xs2 = sv[40];
    real Jrel_np = sv[41];
    real Jrel_p = sv[42];

    #include "ToRORd_fkatp_mixed_endo_mid_epi.common.c"
}

void RHS_RL_cpu(real *a_, real *b_, const real *sv, real *rDY_, real stim_current, real dt, real mapping) {

    // Get the celltype for the current cell
    real celltype = mapping;

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real v = sv[0];
    real CaMKt = sv[1];
    real cass = sv[2];
    real nai = sv[3];
    real nass = sv[4];
    real ki = sv[5];
    real kss = sv[6];
    real cansr = sv[7];
    real cajsr = sv[8];
    real cai = sv[9];
    real m = sv[10];
    real h = sv[11];
    real j = sv[12];
    real hp = sv[13];
    real jp = sv[14];
    real mL = sv[15];
    real hL = sv[16];
    real hLp = sv[17];
    real a = sv[18];
    real iF = sv[19];
    real iS = sv[20];
    real ap = sv[21];
    real iFp = sv[22];
    real iSp = sv[23];
    real d = sv[24];
    real ff = sv[25];
    real fs = sv[26];
    real fcaf = sv[27];
    real fcas = sv[28];
    real jca = sv[29];
    real ffp = sv[30];
    real fcafp = sv[31];
    real nca = sv[32];
    real nca_i = sv[33];
    real ikr_c0 = sv[34];
    real ikr_c1 = sv[35];
    real ikr_c2 = sv[36];
    real ikr_i = sv[37];
    real ikr_o = sv[38];
    real xs1 = sv[39];
    real xs2 = sv[40];
    real Jrel_np = sv[41];
    real Jrel_p = sv[42];

    #include "ToRORd_fkatp_mixed_endo_mid_epi_RL.common.c"
}
