#include "trovato_2019.h"
#include <stdlib.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using Trovato_2019 CPU model\n");

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

        // TODO: Calculate the steady-state for 100 pulses beats (BCL=1000ms)
        // Steady-state 40 pulses (BCL=1000ms)
        sv[0] = -86.7099;
        sv[1] = 0.005431;
        sv[2] = 0.000104;
        sv[3] = 8.25533;
        sv[4] = 8.25502;
        sv[5] = 8.25503;
        sv[6] = 143.743;
        sv[7] = 143.744;
        sv[8] = 143.744;
        sv[9] = 4.4e-05;
        sv[10] = 0.000103;
        sv[11] = 1.26947;
        sv[12] = 1.25254;
        sv[13] = 1.27103;
        sv[14] = 1.1e-05;
        sv[15] = 0;
        sv[16] = 0.006303;
        sv[17] = 0.789469;
        sv[18] = 0.789392;
        sv[19] = 0.791301;
        sv[20] = 0.580955;
        sv[21] = 0.791719;
        sv[22] = 0.000241;
        sv[23] = 0.463851;
        sv[24] = 0.239936;
        sv[25] = 0.000272;
        sv[26] = 0.646362;
        sv[27] = 0.98999;
        sv[28] = 0;
        sv[29] = 1;
        sv[30] = 0.926919;
        sv[31] = 1;
        sv[32] = 1;
        sv[33] = 0.999976;
        sv[34] = 1;
        sv[35] = 1;
        sv[36] = 0.005885;
        sv[37] = 0.000303;
        sv[38] = 0.994251;
        sv[39] = 0.000367;
        sv[40] = 0.566131;
        sv[41] = 0.189842;
        sv[42] = 0.000222;
        sv[43] = 0.233515;
        sv[44] = 0.997077;
        sv[45] = 0.471259;
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

#pragma omp parallel for private(sv_id)
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            //solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver);
            //solve_forward_euler_cpu_adpt_2(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver);
            solve_rush_larsen_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver);
        }
        else {
            for (int j = 0; j < num_steps; ++j) {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current) {

    // -------------------------------------------------------------------------------------------
    // MODEL SPECIFIC:
    // set the variables which are non-linear and hodkin-huxley type
    const real TOLERANCE = 1e-8;
    bool is_rush_larsen[NEQ];
    for (int i = 0; i < NEQ; i++) {
        is_rush_larsen[i] = ((i >= 16 && i <= 35) || (i >= 37 && i <= 44)) ? true : false;        
    }
    // -------------------------------------------------------------------------------------------
    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    // Compute 'a', 'b' coefficients alongside 'rhs'
    real a[NEQ], b[NEQ];
    RHS_RL_cpu(a, b, sv, rDY, stim_current, dt);

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

// Sachetto`s version
void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver) {

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

    RHS_cpu(sv, rDY, stim_curr, *dt);
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
        RHS_cpu(sv, rDY, stim_curr, *dt);
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

// Jhonny`s version
void solve_forward_euler_cpu_adpt_2(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver) {
    const real _beta_safety_ = 0.8;
    const real TOLERANCE = 1e-8;
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

    RHS_cpu(sv, rDY, stim_curr, *dt);
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
        RHS_cpu(sv, rDY, stim_curr, *dt);
        *time_new -= *dt; // step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < numEDO; i++) {
            _k2__[i] = rDY[i];
            real f = (_k1__[i] + _k2__[i]) * 0.5;
            real y_2nd_order = edos_old_aux_[i] + (*dt) * f;
            auxError = (fabs(y_2nd_order) < TOLERANCE) ? fabs(edos_new_euler_[i] - TOLERANCE) : fabs( (y_2nd_order - edos_new_euler_[i])/(y_2nd_order) );
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

void solve_rush_larsen_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver) {
    
    // TODO: Remove this boolean array and write the full algebraics instead ...
    // -------------------------------------------------------------------------------------------
    // MODEL SPECIFIC:
    // set the variables which are non-linear and hodkin-huxley type
    const real TOLERANCE = 1e-8;
    bool is_rush_larsen[NEQ];
    for (int i = 0; i < NEQ; i++) {
        is_rush_larsen[i] = ((i >= 16 && i <= 35) || (i >= 37 && i <= 44)) ? true : false;        
    }
    // -------------------------------------------------------------------------------------------

    //const real _beta_safety_ = 0.85;
    //const real rel_tol = 1.445;
    int numEDO = NEQ;

    real rDY[numEDO];

    // initializes the variables
    solver->ode_previous_dt[sv_id] = solver->ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k2__ = (real *)malloc(sizeof(real) * numEDO);
    real *a_ = (real *)malloc(sizeof(real) * numEDO);
    real *a_new = (real *)malloc(sizeof(real) * numEDO);
    real *b_ = (real *)malloc(sizeof(real) * numEDO);
    real *b_new = (real *)malloc(sizeof(real) * numEDO);
    real *_k_aux__, *_a_aux__, *_b_aux__;

    real *dt = &solver->ode_dt[sv_id];
    real *time_new = &solver->ode_time_new[sv_id];
    real *previous_dt = &solver->ode_previous_dt[sv_id];

    // Keep 'dt' inside the adaptive interval
    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    RHS_RL_cpu(a_, b_, sv, rDY, stim_curr, *dt);
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
                edos_new_euler_[i] = (a_[i] < TOLERANCE) ? edos_old_aux_[i] + (edos_old_aux_[i] * a_[i] + b_[i])*(*dt) : \
                                                  exp(a_[i]*(*dt))*(edos_old_aux_[i] + (b_[i] / a_[i])) - (b_[i] / a_[i]);
            else
                edos_new_euler_[i] = _k1__[i] * *dt + edos_old_aux_[i];
            // steps ahead to compute the rk2 method
            sv[i] = edos_new_euler_[i];
        }

        *time_new += *dt;
        RHS_RL_cpu(a_new, b_new, sv, rDY, stim_curr, *dt);
        *time_new -= *dt; // step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < numEDO; i++) {
            _k2__[i] = rDY[i];
            if (is_rush_larsen[i]) {
                real as = (a_[i] + a_new[i]) * 0.5;
                real bs = (b_[i] + b_new[i]) * 0.5;
                real y_2nd_order = (fabs(as) < TOLERANCE) ? edos_old_aux_[i] + (*dt) * (edos_old_aux_[i]*as + bs) : \
                                                       exp(as*(*dt))*(edos_old_aux_[i] + (bs/as)) - (bs/as);
                auxError = (fabs(y_2nd_order) < TOLERANCE) ? fabs(edos_new_euler_[i] - TOLERANCE) : \
                                                        fabs( (y_2nd_order - edos_new_euler_[i])/(y_2nd_order) );
                greatestError = (auxError > greatestError) ? auxError : greatestError;
            }
            else {
                real f = (_k1__[i] + _k2__[i]) * 0.5;
                real y_2nd_order = edos_old_aux_[i] + (*dt) * f;
                auxError = (fabs(y_2nd_order) < TOLERANCE) ? fabs(edos_new_euler_[i] - TOLERANCE) : \
                                                        fabs( (y_2nd_order - edos_new_euler_[i])/(y_2nd_order) );
                greatestError = (auxError > greatestError) ? auxError : greatestError;
            }
        }
        /// adapt the time step
        greatestError += __tiny_;
        *previous_dt = *dt;
        /// adapt the time step
        //*dt = _beta_safety_ * (*dt) * sqrt(1.0f / greatestError);   // Sachetto`s formula
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

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    const real V = sv[0];
    const real CaMKt = sv[1];
    const real cass = sv[2];
    const real nai = sv[3];
    const real nasl = sv[4];
    const real nass = sv[5];
    const real ki = sv[6];
    const real kss = sv[7];
    const real ksl = sv[8];
    const real cai = sv[9];
    const real casl = sv[10];
    const real cansr = sv[11];
    const real cajsr = sv[12];
    const real cacsr = sv[13];
    const real Jrel1 = sv[14];
    const real Jrel2 = sv[15];
    const real m = sv[16];
    const real hf = sv[17];
    const real hs = sv[18];
    const real j = sv[19];
    const real hsp = sv[20];
    const real jp = sv[21];
    const real mL = sv[22];
    const real hL = sv[23];
    const real hLp = sv[24];
    const real a = sv[25];
    const real i1 = sv[26];
    const real i2 = sv[27];
    const real d = sv[28];
    const real ff = sv[29];
    const real fs = sv[30];
    const real fcaf = sv[31];
    const real fcas = sv[32];
    const real jca = sv[33];
    const real ffp = sv[34];
    const real fcafp = sv[35];
    const real nca = sv[36];
    const real b = sv[37];
    const real g = sv[38];
    const real xrf = sv[39];
    const real xrs = sv[40];
    const real xs1 = sv[41];
    const real xs2 = sv[42];
    const real y = sv[43];
    const real xk1 = sv[44];
    const real u = sv[45];

    #include "trovato_2019_common.inc.c"
}

void RHS_RL_cpu (real *a_, real *b_, real *sv, real *rDY_, real stim_current, real dt) {
    
    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    const real V = sv[0];
    const real CaMKt = sv[1];
    const real cass = sv[2];
    const real nai = sv[3];
    const real nasl = sv[4];
    const real nass = sv[5];
    const real ki = sv[6];
    const real kss = sv[7];
    const real ksl = sv[8];
    const real cai = sv[9];
    const real casl = sv[10];
    const real cansr = sv[11];
    const real cajsr = sv[12];
    const real cacsr = sv[13];
    const real Jrel1 = sv[14];
    const real Jrel2 = sv[15];
    const real m = sv[16];
    const real hf = sv[17];
    const real hs = sv[18];
    const real j = sv[19];
    const real hsp = sv[20];
    const real jp = sv[21];
    const real mL = sv[22];
    const real hL = sv[23];
    const real hLp = sv[24];
    const real a = sv[25];
    const real i1 = sv[26];
    const real i2 = sv[27];
    const real d = sv[28];
    const real ff = sv[29];
    const real fs = sv[30];
    const real fcaf = sv[31];
    const real fcas = sv[32];
    const real jca = sv[33];
    const real ffp = sv[34];
    const real fcafp = sv[35];
    const real nca = sv[36];
    const real b = sv[37];
    const real g = sv[38];
    const real xrf = sv[39];
    const real xrs = sv[40];
    const real xs1 = sv[41];
    const real xs2 = sv[42];
    const real y = sv[43];
    const real xk1 = sv[44];
    const real u = sv[45];

    #include "trovato_2019_RL_common.inc.c"
}
