#include "ToRORd_Land_mixed_endo_mid_epi.h"
#include <stdlib.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using ToRORd_Land_2020 CPU model\n");

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
        struct extra_data_for_torord_land *extra_data = (struct extra_data_for_torord_land*)solver->ode_extra_data;
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

            // Default initial conditions for ENDO cell (from original Matlab script)
            sv[0] = -8.863699e+01;
            sv[1] = 1.189734e+01;
            sv[2] = 1.189766e+01;
            sv[3] = 1.412345e+02;
            sv[4] = 1.412344e+02;
            sv[5] = 7.267473e-05;
            sv[6] = 6.337870e-05;
            sv[7] = 1.532653e+00;
            sv[8] = 1.533946e+00;
            sv[9] = 8.280078e-04;
            sv[10] = 6.665272e-01;
            sv[11] = 8.260208e-01;
            sv[12] = 8.260560e-01;
            sv[13] = 8.258509e-01;
            sv[14] = 1.668686e-04;
            sv[15] = 5.228306e-01;
            sv[16] = 2.859696e-01;
            sv[17] = 9.591370e-04;
            sv[18] = 9.996012e-01;
            sv[19] = 5.934016e-01;
            sv[20] = 4.886961e-04;
            sv[21] = 9.996011e-01;
            sv[22] = 6.546687e-01;
            sv[23] = 9.500075e-32;
            sv[24] = 1.000000e+00;
            sv[25] = 9.392580e-01;
            sv[26] = 1.000000e+00;
            sv[27] = 9.998984e-01;
            sv[28] = 9.999783e-01;
            sv[29] = 4.448162e-04;
            sv[30] = 7.550725e-04;
            sv[31] = 1.000000e+00;
            sv[32] = 1.000000e+00;
            sv[33] = 2.424047e-01;
            sv[34] = 1.795377e-04;
            sv[35] = -6.883086e-25;
            sv[36] = 1.117498e-02;
            sv[37] = 9.980366e-01;
            sv[38] = 8.588018e-04;
            sv[39] = 7.097447e-04;
            sv[40] = 3.812617e-04;
            sv[41] = 1.357116e-05;
            sv[42] = 2.302525e-23;
            sv[43] = 1.561941e-04;
            sv[44] = 2.351289e-04;
            sv[45] = 8.077631e-03;
            sv[46] = 9.993734e-01;
            sv[47] = 0.000000e+00;
            sv[48] = 0.000000e+00;
	        sv[49] = 0.000000e+00;

            // Default initial conditions for MID cell (from original Matlab script)
            //sv[0] = -8.953800e+01;
            //sv[1] = 1.492920e+01;
            //sv[2] = 1.492967e+01;
            //sv[3] = 1.448447e+02;
            //sv[4] = 1.448447e+02;
            //sv[5] = 7.502288e-05;
            //sv[6] = 6.107636e-05;
            //sv[7] = 1.790435e+00;
            //sv[8] = 1.794842e+00;
            //sv[9] = 6.819365e-04;
            //sv[10] = 6.953807e-01;
            //sv[11] = 8.434888e-01;
            //sv[12] = 8.435208e-01;
            //sv[13] = 8.432262e-01;
            //sv[14] = 1.406211e-04;
            //sv[15] = 5.453149e-01;
            //sv[16] = 2.924967e-01;
            //sv[17] = 9.026127e-04;
            //sv[18] = 9.996593e-01;
            //sv[19] = 5.631197e-01;
            //sv[20] = 4.598833e-04;
            //sv[21] = 9.996593e-01;
            //sv[22] = 6.236964e-01;
            //sv[23] = -1.314189e-33;
            //sv[24] = 1.000000e+00;
            //sv[25] = 9.204086e-01;
            //sv[26] = 1.000000e+00;
            //sv[27] = 9.997620e-01;
            //sv[28] = 9.999625e-01;
            //sv[29] = 3.853595e-04;
            //sv[30] = 8.535292e-04;
            //sv[31] = 1.000000e+00;
            //sv[32] = 1.000000e+00;
            //sv[33] = 2.664151e-01;
            //sv[34] = 1.623107e-04;
            //sv[35] = 1.209762e-24;
            //sv[36] = 1.782437e-02;
            //sv[37] = 9.979720e-01;
            //sv[38] = 8.053991e-04;
            //sv[39] = 6.781800e-04;
            //sv[40] = 5.265363e-04;
            //sv[41] = 1.789565e-05;
            //sv[42] = 7.059162e-23;
            //sv[43] = 1.670654e-04;
            //sv[44] = 2.506794e-04;
            //sv[45] = 8.602625e-03;
            //sv[46] = 9.993314e-01;
            //sv[47] = 0.000000e+00;
            //sv[48] = 0.000000e+00;

            // Default initial conditions for EPI cell (from original Matlab script)
            //sv[0] = -8.904628e+01;
            //sv[1] = 1.272190e+01;
            //sv[2] = 1.272220e+01;
            //sv[3] = 1.422490e+02;
            //sv[4] = 1.422489e+02;
            //sv[5] = 6.541058e-05;
            //sv[6] = 5.684431e-05;
            //sv[7] = 1.809117e+00;
            //sv[8] = 1.809702e+00;
            //sv[9] = 7.581821e-04;
            //sv[10] = 6.798398e-01;
            //sv[11] = 8.341502e-01;
            //sv[12] = 8.341883e-01;
            //sv[13] = 8.340817e-01;
            //sv[14] = 1.543877e-04;
            //sv[15] = 5.382951e-01;
            //sv[16] = 3.027694e-01;
            //sv[17] = 9.330351e-04;
            //sv[18] = 9.996287e-01;
            //sv[19] = 9.996262e-01;
            //sv[20] = 4.753907e-04;
            //sv[21] = 9.996287e-01;
            //sv[22] = 9.996285e-01;
            //sv[23] = 1.742134e-37;
            //sv[24] = 1.000000e+00;
            //sv[25] = 9.479522e-01;
            //sv[26] = 1.000000e+00;
            //sv[27] = 9.999327e-01;
            //sv[28] = 9.999829e-01;
            //sv[29] = 2.915447e-04;
            //sv[30] = 5.026045e-04;
            //sv[31] = 1.000000e+00;
            //sv[32] = 1.000000e+00;
            //sv[33] = 2.288155e-01;
            //sv[34] = 1.714978e-04;
            //sv[35] = -1.131190e-26;
            //sv[36] = 1.295052e-02;
            //sv[37] = 9.981944e-01;
            //sv[38] = 8.342321e-04;
            //sv[39] = 6.838658e-04;
            //sv[40] = 2.778785e-04;
            //sv[41] = 9.667759e-06;
            //sv[42] = 8.169304e-24;
            //sv[43] = 1.259996e-04;
            //sv[44] = 1.899522e-04;
            //sv[45] = 6.551494e-03;
            //sv[46] = 9.994940e-01;
            //sv[47] = 0.000000e+00;
            //sv[48] = 0.000000e+00;
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
    int num_extra_parameters = 20;
    real extra_par[num_extra_parameters];
    real *transmurality = NULL;
    if (ode_solver->ode_extra_data) {
        struct extra_data_for_torord_land *extra_data = (struct extra_data_for_torord_land*)ode_solver->ode_extra_data;
        extra_par[0]  = extra_data->INa_Multiplier; 
        extra_par[1]  = extra_data->INaL_Multiplier;
        extra_par[2]  = extra_data->INaCa_Multiplier;
        extra_par[3]  = extra_data->INaK_Multiplier;
        extra_par[4]  = extra_data->INab_Multiplier; 
        extra_par[5]  = extra_data->Ito_Multiplier;
        extra_par[6]  = extra_data->IKr_Multiplier; 
        extra_par[7]  = extra_data->IKs_Multiplier; 
        extra_par[8]  = extra_data->IK1_Multiplier;
        extra_par[9]  = extra_data->IKb_Multiplier;
        extra_par[10]  = extra_data->IKCa_Multiplier;
        extra_par[11] = extra_data->ICaL_Multiplier;  
        extra_par[12] = extra_data->ICab_Multiplier;  
        extra_par[13] = extra_data->IpCa_Multiplier;
        extra_par[14] = extra_data->ICaCl_Multiplier; 
        extra_par[15] = extra_data->IClb_Multiplier;
        extra_par[16] = extra_data->Jrel_Multiplier;
        extra_par[17] = extra_data->Jup_Multiplier;
        extra_par[18] = extra_data->aCaMK_Multiplier;
        extra_par[19] = extra_data->taurelp_Multiplier;
        transmurality = extra_data->transmurality;
    }
    else {
        // Default: initialize all current modifiers
        for (uint32_t i = 0; i < num_extra_parameters; i++) {
            if (i == 10)
                extra_par[i] = 0.0;
            else
                extra_par[i] = 1.0;
        }
    }

    // Solve the ODEs
    OMP(parallel for private(sv_id))
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            if (ode_solver->ode_extra_data) {
                solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], transmurality[i], current_t + dt, sv_id, ode_solver, extra_par);
            }
            else {
                solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], 0.0, current_t + dt, sv_id, ode_solver, extra_par);
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
    SOLVE_EQUATION_EULER_CPU(1);        // nai    
    SOLVE_EQUATION_EULER_CPU(2);        // nass 
    SOLVE_EQUATION_EULER_CPU(3);        // ki 
    SOLVE_EQUATION_EULER_CPU(4);        // kss 
    SOLVE_EQUATION_EULER_CPU(5);        // cai   
    SOLVE_EQUATION_EULER_CPU(6);        // cass  
    SOLVE_EQUATION_EULER_CPU(7);        // cansr
    SOLVE_EQUATION_EULER_CPU(8);        // cajsr
    SOLVE_EQUATION_RUSH_LARSEN_CPU(9);  // m
    SOLVE_EQUATION_RUSH_LARSEN_CPU(10); // hp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(11); // h
    SOLVE_EQUATION_RUSH_LARSEN_CPU(12); // j
    SOLVE_EQUATION_RUSH_LARSEN_CPU(13); // jp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(14); // mL
    SOLVE_EQUATION_RUSH_LARSEN_CPU(15); // hL
    SOLVE_EQUATION_RUSH_LARSEN_CPU(16); // hLp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(17); // a
    SOLVE_EQUATION_RUSH_LARSEN_CPU(18); // iF
    SOLVE_EQUATION_RUSH_LARSEN_CPU(19); // iS
    SOLVE_EQUATION_RUSH_LARSEN_CPU(20); // ap
    SOLVE_EQUATION_RUSH_LARSEN_CPU(21); // iFp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(22); // iSp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(23); // d
    SOLVE_EQUATION_RUSH_LARSEN_CPU(24); // ff
    SOLVE_EQUATION_RUSH_LARSEN_CPU(25); // fs
    SOLVE_EQUATION_RUSH_LARSEN_CPU(26); // fcaf
    SOLVE_EQUATION_RUSH_LARSEN_CPU(27); // fcas
    SOLVE_EQUATION_RUSH_LARSEN_CPU(28); // jca
    SOLVE_EQUATION_EULER_CPU(29);       // nca
    SOLVE_EQUATION_EULER_CPU(30);       // nca_i
    SOLVE_EQUATION_RUSH_LARSEN_CPU(31); // ffp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(32); // fcafp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(33); // xs1
    SOLVE_EQUATION_RUSH_LARSEN_CPU(34); // xs2
    SOLVE_EQUATION_RUSH_LARSEN_CPU(35); // Jrel_np
    SOLVE_EQUATION_EULER_CPU(36);       // CaMKt
    SOLVE_EQUATION_EULER_CPU(37);       // ikr_c0
    SOLVE_EQUATION_EULER_CPU(38);       // ikr_c1
    SOLVE_EQUATION_EULER_CPU(39);       // ikr_c2
    SOLVE_EQUATION_EULER_CPU(40);       // ikr_o
    SOLVE_EQUATION_EULER_CPU(41);       // ikr_i
    SOLVE_EQUATION_RUSH_LARSEN_CPU(42); // Jrel_p
    // ---------------------------------------------------
    // Land-Niederer
    SOLVE_EQUATION_EULER_CPU(43);       // XS
    SOLVE_EQUATION_EULER_CPU(44);       // XW
    SOLVE_EQUATION_EULER_CPU(45);       // Ca_TRPN
    SOLVE_EQUATION_EULER_CPU(46);       // TmBlocked
    SOLVE_EQUATION_EULER_CPU(47);       // ZETAS
    SOLVE_EQUATION_EULER_CPU(48);       // ZETAW 
    SOLVE_EQUATION_CONSTANT_CPU(49);	// Ta
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

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real const *extra_params) {

    // Current modifiers
    real INa_Multiplier = extra_params[0];   
    real INaL_Multiplier = extra_params[1];  
    real INaCa_Multiplier = extra_params[2];  
    real INaK_Multiplier = extra_params[3];  
    real INab_Multiplier = extra_params[4];   
    real Ito_Multiplier = extra_params[5];  
    real IKr_Multiplier = extra_params[6];   
    real IKs_Multiplier = extra_params[7];   
    real IK1_Multiplier = extra_params[8];  
    real IKb_Multiplier = extra_params[9];  
    real IKCa_Multiplier = extra_params[10];  
    real ICaL_Multiplier = extra_params[11];   
    real ICab_Multiplier = extra_params[12];   
    real IpCa_Multiplier = extra_params[13]; 
    real ICaCl_Multiplier = extra_params[14];  
    real IClb_Multiplier = extra_params[15]; 
    real Jrel_Multiplier = extra_params[16]; 
    real Jup_Multiplier = extra_params[17]; 
    real aCaMK_Multiplier = extra_params[18]; 
    real taurelp_Multiplier = extra_params[19]; 

    // Get the celltype for the current cell
    real celltype = transmurality;

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables (same order as the original Matlab script)
    real v = sv[0];
    real nai = sv[1];
    real nass = sv[2];
    real ki = sv[3];
    real kss = sv[4];
    real cai = sv[5];
    real cass = sv[6];
    real cansr = sv[7];
    real cajsr = sv[8];
    real m = sv[9];
    real hp = sv[10];
    real h = sv[11];
    real j = sv[12];
    real jp = sv[13];
    real mL = sv[14];
    real hL = sv[15];
    real hLp = sv[16];
    real a = sv[17];
    real iF = sv[18];
    real iS = sv[19];
    real ap = sv[20];
    real iFp = sv[21];
    real iSp = sv[22];

    // ical
    real d = sv[23];
    real ff = sv[24];
    real fs = sv[25];
    real fcaf = sv[26];
    real fcas = sv[27];
    real jca = sv[28];
    real nca = sv[29];
    real nca_i = sv[30];
    real ffp = sv[31];
    real fcafp = sv[32];

    real xs1 = sv[33];
    real xs2 = sv[34];
    real Jrel_np = sv[35];
    real CaMKt = sv[36];

    // new MM ICaL states
    real ikr_c0 = sv[37];
    real ikr_c1 = sv[38];
    real ikr_c2 = sv[39];
    real ikr_o = sv[40];
    real ikr_i = sv[41];
    real Jrel_p = sv[42];

    const real cli = 24;   // Intracellular Cl  [mM]
    const real clo = 150;  // Extracellular Cl  [mM]
// -----------------------------------------------------
    // Land-Niederer model
    real XS = fmaxf(0,sv[43]);
    real XW = fmaxf(0,sv[44]);
    real Ca_TRPN = fmaxf(0,sv[45]);
    real TmBlocked = sv[46];
    real ZETAS = sv[47];
    real ZETAW = sv[48];
    real TA = sv[49];

    #include "ToRORd_Land_mixed_endo_mid_epi.common.c"
}

void RHS_RL_cpu(real *a_, real *b_, const real *sv, real *rDY_, real stim_current, real dt, \
                real transmurality, real const *extra_params) {

    // Current modifiers
    real INa_Multiplier = extra_params[0];
    real INaL_Multiplier = extra_params[1];  
    real INaCa_Multiplier = extra_params[2];  
    real INaK_Multiplier = extra_params[3];  
    real INab_Multiplier = extra_params[4];   
    real Ito_Multiplier = extra_params[5];  
    real IKr_Multiplier = extra_params[6];   
    real IKs_Multiplier = extra_params[7];   
    real IK1_Multiplier = extra_params[8];  
    real IKb_Multiplier = extra_params[9];  
    real IKCa_Multiplier = extra_params[10];  
    real ICaL_Multiplier = extra_params[11];   
    real ICab_Multiplier = extra_params[12];   
    real IpCa_Multiplier = extra_params[13]; 
    real ICaCl_Multiplier = extra_params[14];  
    real IClb_Multiplier = extra_params[15]; 
    real Jrel_Multiplier = extra_params[16]; 
    real Jup_Multiplier = extra_params[17]; 
    real aCaMK_Multiplier = extra_params[18]; 
    real taurelp_Multiplier = extra_params[19]; 

    // Get the celltype for the current cell
    real celltype = transmurality;

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables (same order as the original Matlab script)
    real v = sv[0];
    real nai = sv[1];
    real nass = sv[2];
    real ki = sv[3];
    real kss = sv[4];
    real cai = sv[5];
    real cass = sv[6];
    real cansr = sv[7];
    real cajsr = sv[8];
    real m = sv[9];
    real hp = sv[10];
    real h = sv[11];
    real j = sv[12];
    real jp = sv[13];
    real mL = sv[14];
    real hL = sv[15];
    real hLp = sv[16];
    real a = sv[17];
    real iF = sv[18];
    real iS = sv[19];
    real ap = sv[20];
    real iFp = sv[21];
    real iSp = sv[22];

    // ical
    real d = sv[23];
    real ff = sv[24];
    real fs = sv[25];
    real fcaf = sv[26];
    real fcas = sv[27];
    real jca = sv[28];
    real nca = sv[29];
    real nca_i = sv[30];
    real ffp = sv[31];
    real fcafp = sv[32];

    real xs1 = sv[33];
    real xs2 = sv[34];
    real Jrel_np = sv[35];
    real CaMKt = sv[36];

    // new MM ICaL states
    real ikr_c0 = sv[37];
    real ikr_c1 = sv[38];
    real ikr_c2 = sv[39];
    real ikr_o = sv[40];
    real ikr_i = sv[41];
    real Jrel_p = sv[42];

    const real cli = 24;   // Intracellular Cl  [mM]
    const real clo = 150;  // Extracellular Cl  [mM]
// -----------------------------------------------------
    real XS = fmaxf(0,sv[43]);
    real XW = fmaxf(0,sv[44]);
    real Ca_TRPN = fmaxf(0,sv[45]);
    real TmBlocked = sv[46];
    real ZETAS = sv[47];
    real ZETAW = sv[48];

    real TA = sv[49];

    #include "ToRORd_Land_mixed_endo_mid_epi_RL.common.c"
}
