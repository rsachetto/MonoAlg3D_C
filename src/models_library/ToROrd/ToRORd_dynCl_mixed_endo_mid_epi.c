#include "ToRORd_dynCl_mixed_endo_mid_epi.h"
#include <stdlib.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using ToRORd_dynCl_2020 CPU model\n");

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
    real *transmurality = NULL;
    if(solver->ode_extra_data) {
        struct extra_data_for_torord *extra_data = (struct extra_data_for_torord*)solver->ode_extra_data;
        initial_endo = extra_data->initial_ss_endo;
        initial_epi = extra_data->initial_ss_epi;
        initial_mid = extra_data->initial_ss_mid;
        transmurality = extra_data->transmurality;

        OMP(parallel for)
        for(uint32_t i = 0; i < num_cells; i++){
            
            real *sv = &solver->sv[i * NEQ];

            for (int j = 0; j < NEQ; j++) {
                if (transmurality[i] == ENDO)
                    sv[j] = initial_endo[j];
                else if (transmurality[i] == EPI)
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

            // Initial conditions 200 beats (endocardium cell)
            sv[0] = -9.035192e+01;
            sv[1] = 1.162900e-02;
            sv[2] = 6.500000e-05;
            sv[3] = 1.239893e+01;
            sv[4] = 1.239926e+01;
            sv[5] = 1.482415e+02;
            sv[6] = 1.482414e+02;
            sv[7] = 1.527292e+00;
            sv[8] = 1.524395e+00;
            sv[9] = 7.400000e-05;
            sv[10] = 5.720000e-04;
            sv[11] = 8.579420e-01;
            sv[12] = 8.577990e-01;
            sv[13] = 7.199660e-01;
            sv[14] = 8.575760e-01;
            sv[15] = 1.200000e-04;
            sv[16] = 5.748970e-01;
            sv[17] = 3.250180e-01;
            sv[18] = 8.540000e-04;
            sv[19] = 9.997050e-01;
            sv[20] = 5.959350e-01;
            sv[21] = 4.350000e-04;
            sv[22] = 9.997050e-01;
            sv[23] = 6.589890e-01;
            sv[24] = 0.000000e+00;
            sv[25] = 1.000000e+00;
            sv[26] = 9.343710e-01;
            sv[27] = 1.000000e+00;
            sv[28] = 9.998810e-01;
            sv[29] = 9.999820e-01;
            sv[30] = 1.000000e+00;
            sv[31] = 1.000000e+00;
            sv[32] = 4.830000e-04;
            sv[33] = 8.180000e-04;
            sv[34] = 9.983340e-01;
            sv[35] = 7.600000e-04;
            sv[36] = 6.260000e-04;
            sv[37] = 9.000000e-06;
            sv[38] = 2.720000e-04;
            sv[39] = 2.568150e-01;
            sv[40] = 1.480000e-04;
            sv[41] = 0.000000e+00;
            sv[42] = 0.000000e+00;
            sv[43] = 2.978204e+01;
            sv[44] = 2.978201e+01;

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

    // Get the extra parameters
    int num_extra_parameters = 17;
    real extra_par[num_extra_parameters];
    real *transmurality = NULL;
    if (ode_solver->ode_extra_data) {
        struct extra_data_for_torord *extra_data = (struct extra_data_for_torord*)ode_solver->ode_extra_data;
        extra_par[0]  = extra_data->INa_Multiplier; 
        extra_par[1]  = extra_data->ICaL_Multiplier;
        extra_par[2]  = extra_data->Ito_Multiplier;
        extra_par[3]  = extra_data->INaL_Multiplier;
        extra_par[4]  = extra_data->IKr_Multiplier; 
        extra_par[5]  = extra_data->IKs_Multiplier; 
        extra_par[6]  = extra_data->IK1_Multiplier; 
        extra_par[7]  = extra_data->IKb_Multiplier; 
        extra_par[8]  = extra_data->INaCa_Multiplier;
        extra_par[9]  = extra_data->INaK_Multiplier;  
        extra_par[9]  = extra_data->INab_Multiplier;  
        extra_par[10] = extra_data->ICab_Multiplier;  
        extra_par[11] = extra_data->IpCa_Multiplier;  
        extra_par[12] = extra_data->ICaCl_Multiplier;
        extra_par[13] = extra_data->IClb_Multiplier; 
        extra_par[15] = extra_data->Jrel_Multiplier; 
        extra_par[16] = extra_data->Jup_Multiplier;
        transmurality = extra_data->transmurality;
    }
    else {
        extra_par[0]  = 1.0; 
        extra_par[1]  = 1.0;
        extra_par[2]  = 1.0;
        extra_par[3]  = 1.0;
        extra_par[4]  = 1.0;
        extra_par[5]  = 1.0;
        extra_par[6]  = 1.0; 
        extra_par[7]  = 1.0; 
        extra_par[8]  = 1.0;
        extra_par[9]  = 1.0;
        extra_par[9]  = 1.0; 
        extra_par[10] = 1.0;  
        extra_par[11] = 1.0; 
        extra_par[12] = 1.0;
        extra_par[13] = 1.0;
        extra_par[15] = 1.0;
        extra_par[16] = 1.0;
    }

    OMP(parallel for private(sv_id))
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            if (ode_solver->ode_extra_data) {
                //solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], transmurality[i], current_t + dt, sv_id, ode_solver, extra_par);
                solve_rush_larsen_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], transmurality[i], current_t + dt, sv_id, ode_solver, extra_par);
            }
            else {
                //solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], 0.0, current_t + dt, sv_id, ode_solver, extra_par);
                solve_rush_larsen_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], 0.0, current_t + dt, sv_id, ode_solver, extra_par);
            }
        }
        else {
            for (int j = 0; j < num_steps; ++j) {
                if (ode_solver->ode_extra_data) {
                    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], transmurality[i], extra_par);
                }
                else {
                    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], 0.0, extra_par);
                }
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current, real transmurality, real const *extra_params) {

    const real TOLERANCE = 1e-8;
    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    // Compute 'a', 'b' coefficients alongside 'rhs'
    real a[NEQ], b[NEQ];
    RHS_RL_cpu(a, b, sv, rDY, stim_current, dt, transmurality, extra_params);

    // Solve variables based on its type:
    //  Non-linear = Euler
    //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
    SOLVE_EQUATION_EULER_CPU(0);        // v        
    SOLVE_EQUATION_EULER_CPU(1);        // CaMKt    
    SOLVE_EQUATION_EULER_CPU(2);        // cass 
    SOLVE_EQUATION_EULER_CPU(3);        // nai  
    SOLVE_EQUATION_EULER_CPU(4);        // nass 
    SOLVE_EQUATION_EULER_CPU(5);        // ki   
    SOLVE_EQUATION_EULER_CPU(6);        // kss  
    SOLVE_EQUATION_EULER_CPU(7);        // cansr
    SOLVE_EQUATION_EULER_CPU(8);        // cajsr
    SOLVE_EQUATION_EULER_CPU(9);        // cai
    SOLVE_EQUATION_RUSH_LARSEN_CPU(10); // m
    SOLVE_EQUATION_RUSH_LARSEN_CPU(11); // h
    SOLVE_EQUATION_RUSH_LARSEN_CPU(12); // j
    SOLVE_EQUATION_RUSH_LARSEN_CPU(13); // hp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(14); // jp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(15); // mL
    SOLVE_EQUATION_RUSH_LARSEN_CPU(16); // hL
    SOLVE_EQUATION_RUSH_LARSEN_CPU(17); // hLp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(18); // a
    SOLVE_EQUATION_RUSH_LARSEN_CPU(19); // iF
    SOLVE_EQUATION_RUSH_LARSEN_CPU(20); // iS
    SOLVE_EQUATION_RUSH_LARSEN_CPU(21); // ap
    SOLVE_EQUATION_RUSH_LARSEN_CPU(22); // iFp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(23); // iSp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(24); // d
    SOLVE_EQUATION_RUSH_LARSEN_CPU(25); // ff
    SOLVE_EQUATION_RUSH_LARSEN_CPU(26); // fs
    SOLVE_EQUATION_RUSH_LARSEN_CPU(27); // fcaf
    SOLVE_EQUATION_RUSH_LARSEN_CPU(28); // fcas
    SOLVE_EQUATION_RUSH_LARSEN_CPU(29); // jca
    SOLVE_EQUATION_RUSH_LARSEN_CPU(30); // ffp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(31); // fcafp
    SOLVE_EQUATION_EULER_CPU(32);       // nca
    SOLVE_EQUATION_EULER_CPU(33);       // nca_i
    SOLVE_EQUATION_EULER_CPU(34);       // ikr_c0
    SOLVE_EQUATION_EULER_CPU(35);       // ikr_c1
    SOLVE_EQUATION_EULER_CPU(36);       // ikr_c2
    SOLVE_EQUATION_EULER_CPU(37);       // ikr_i
    SOLVE_EQUATION_EULER_CPU(38);       // ikr_o
    SOLVE_EQUATION_RUSH_LARSEN_CPU(39); // xs1
    SOLVE_EQUATION_RUSH_LARSEN_CPU(40); // xs2
    SOLVE_EQUATION_RUSH_LARSEN_CPU(41); // Jrel_np
    SOLVE_EQUATION_RUSH_LARSEN_CPU(42); // Jrel_p
    SOLVE_EQUATION_EULER_CPU(43);       // cli
    SOLVE_EQUATION_EULER_CPU(44);       // clss
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real transmurality, real final_time, int sv_id, struct ode_solver *solver, real const *extra_params) {

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

    RHS_cpu(sv, rDY, stim_curr, *dt, transmurality, extra_params);
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
        RHS_cpu(sv, rDY, stim_curr, *dt, transmurality, extra_params);
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

void solve_rush_larsen_cpu_adpt(real *sv, real stim_curr, real transmurality, real final_time, int sv_id, struct ode_solver *solver, real const *extra_params) {
    
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

    RHS_RL_cpu(a_, b_, sv, rDY, stim_curr, *dt, transmurality, extra_params);
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

        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(0);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(1);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(2);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(3);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(4);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(5);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(6);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(7);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(8);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(9);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(10);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(11);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(12);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(13);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(14);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(15);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(16);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(17);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(18);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(19);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(20);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(21);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(22);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(23);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(24);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(25);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(26);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(27);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(28);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(29);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(30);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(31);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(32);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(33);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(34);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(35);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(36);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(37);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(38);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(39);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(40);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(41);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(42);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(43);
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(44);

        *time_new += *dt;
        RHS_RL_cpu(a_new, b_new, sv, rDY, stim_curr, *dt, transmurality, extra_params);
        *time_new -= *dt; // step back

        // Compute errors
        double greatestError = 0.0, auxError = 0.0;
        real as, bs, f, y_2nd_order;
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(0);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(1);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(2);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(3);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(4);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(5);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(6);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(7);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(8);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(9);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(10);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(11);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(12);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(13);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(14);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(15);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(16);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(17);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(18);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(19);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(20);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(21);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(22);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(23);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(24);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(25);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(26);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(27);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(28);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(29);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(30);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(31);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(32);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(33);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(34);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(35);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(36);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(37);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(38);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(39);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(40);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(41);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(42);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(43);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(44);
    
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

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real const *extra_params) {

    // Current modifiers
    real INa_Multiplier   = extra_params[0]; 
    real ICaL_Multiplier  = extra_params[1];
    real Ito_Multiplier   = extra_params[2];
    real INaL_Multiplier  = extra_params[3];
    real IKr_Multiplier   = extra_params[4]; 
    real IKs_Multiplier   = extra_params[5]; 
    real IK1_Multiplier   = extra_params[6]; 
    real IKb_Multiplier   = extra_params[7]; 
    real INaCa_Multiplier = extra_params[8];
    real INaK_Multiplier  = extra_params[9];  
    real INab_Multiplier  = extra_params[10];  
    real ICab_Multiplier  = extra_params[11];  
    real IpCa_Multiplier  = extra_params[12];  
    real ICaCl_Multiplier = extra_params[13];
    real IClb_Multiplier  = extra_params[14]; 
    real Jrel_Multiplier  = extra_params[15]; 
    real Jup_Multiplier   = extra_params[16];

    // Get the celltype for the current cell
    real celltype = transmurality;

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
    real cli = sv[43];
    real clss = sv[44];

    #include "ToRORd_dynCl_mixed_endo_mid_epi.common.c"
}

void RHS_RL_cpu(real *a_, real *b_, const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real const *extra_params) {

    // Current modifiers
    real INa_Multiplier   = extra_params[0]; 
    real ICaL_Multiplier  = extra_params[1];
    real Ito_Multiplier   = extra_params[2];
    real INaL_Multiplier  = extra_params[3];
    real IKr_Multiplier   = extra_params[4]; 
    real IKs_Multiplier   = extra_params[5]; 
    real IK1_Multiplier   = extra_params[6]; 
    real IKb_Multiplier   = extra_params[7]; 
    real INaCa_Multiplier = extra_params[8];
    real INaK_Multiplier  = extra_params[9];  
    real INab_Multiplier  = extra_params[10];  
    real ICab_Multiplier  = extra_params[11];  
    real IpCa_Multiplier  = extra_params[12];  
    real ICaCl_Multiplier = extra_params[13];
    real IClb_Multiplier  = extra_params[14]; 
    real Jrel_Multiplier  = extra_params[15]; 
    real Jup_Multiplier   = extra_params[16];

    // Get the celltype for the current cell
    real celltype = transmurality;

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
    real cli = sv[43];
    real clss = sv[44];

    #include "ToRORd_dynCl_mixed_endo_mid_epi_RL.common.c"
}
