#include "ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.h"
#include <stdlib.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using ToRORd_fkatp_2019 CPU model with GKs GKr and tjca adjustments\n");

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
	real *sf_Iks = NULL;
    if(solver->ode_extra_data) {
        struct extra_data_for_torord_gksgkrtjca_twave *extra_data = (struct extra_data_for_torord_gksgkrtjca_twave*)solver->ode_extra_data;
        initial_endo = extra_data->initial_ss_endo;
        initial_epi = extra_data->initial_ss_epi;
        initial_mid = extra_data->initial_ss_mid;
        transmurality = extra_data->transmurality;
		sf_Iks = extra_data->sf_IKs;

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
            // Default initial condition from Matlab (endocardium cell)
            sv[0]  = -88.6369922306458;
            sv[1]  = 11.8973412949238;
            sv[2]  = 11.897661047085;
            sv[3]  = 141.234464714982;
            sv[4]  = 141.234423402713;
            sv[5]  = 7.26747296460659e-05;
            sv[6]  = 6.33786975780735e-05;
            sv[7]  = 1.5326530637197;
            sv[8]  = 1.53394579180493;
            sv[9]  = 0.000828007761976018;
            sv[10] = 0.666527193684116;
            sv[11] = 0.826020806005678;
            sv[12] = 0.826055985895856;
            sv[13] = 0.825850881115628;
            sv[14] = 0.000166868626513013;
            sv[15] = 0.522830604669169;
            sv[16] = 0.285969584294187;
            sv[17] = 0.000959137028030184;
            sv[18] = 0.999601150012565;
            sv[19] = 0.5934016398361;
            sv[20] = 0.000488696137242056;
            sv[21] = 0.999601147267179;
            sv[22] = 0.654668660159696;
            sv[23] = 9.50007519781516e-32;
            sv[24] = 0.999999992317577;
            sv[25] = 0.939258048397962;
            sv[26] = 0.999999992317557;
            sv[27] = 0.999898379647465;
            sv[28] = 0.99997825156004;
            sv[29] = 0.000444816183420527;
            sv[30] = 0.000755072490632667;
            sv[31] = 0.999999992318446;
            sv[32] = 0.999999992318445;
            sv[33] = 0.24240468344952;
            sv[34] = 0.000179537726989804;
            sv[35] = -6.88308558109975e-25;
            sv[36] = 0.0111749845355653;
            sv[37] = 0.998036620213316;
            sv[38] = 0.000858801779013532;
            sv[39] = 0.000709744678350176;
            sv[40] = 0.000381261722195702;
            sv[41] = 1.35711566929992e-05;
            sv[42] = 2.30252452954649e-23;
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
	real *sf_Iks = NULL;
    if (ode_solver->ode_extra_data) {
        struct extra_data_for_torord_gksgkrtjca_twave *extra_data = (struct extra_data_for_torord_gksgkrtjca_twave*)ode_solver->ode_extra_data;
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
        sf_Iks = extra_data->sf_IKs;
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

    OMP(parallel for private(sv_id))
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            if (ode_solver->ode_extra_data) {
                solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], transmurality[i], sf_Iks[i], current_t + dt, sv_id, ode_solver, extra_par);
            }
            else {
			    solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], 0.0, 1.0, current_t + dt, sv_id, ode_solver, extra_par);
            }
        }
        else {
            for (int j = 0; j < num_steps; ++j) {
                if (ode_solver->ode_extra_data) {
                    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], transmurality[i], sf_Iks[i], extra_par);
                }
                else {
				    solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], 0.0, 1.0, extra_par);
                }
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current, real transmurality, real sf_Iks, real const *extra_params) {

    const real TOLERANCE = 1e-8;
    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    // Compute 'a', 'b' coefficients alongside 'rhs'
    real a[NEQ], b[NEQ];
    RHS_RL_cpu(a, b, sv, rDY, stim_current, dt, transmurality, sf_Iks, extra_params);

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
    SOLVE_EQUATION_EULER_CPU(36);        // CaMKt 
    SOLVE_EQUATION_EULER_CPU(37);       // ikr_c0
    SOLVE_EQUATION_EULER_CPU(38);       // ikr_c1
    SOLVE_EQUATION_EULER_CPU(39);       // ikr_c2
    SOLVE_EQUATION_EULER_CPU(40);       // ikr_o
    SOLVE_EQUATION_EULER_CPU(41);       // ikr_i
    SOLVE_EQUATION_RUSH_LARSEN_CPU(42); // Jrel_p
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real transmurality, real sf_Iks, real final_time, int sv_id, struct ode_solver *solver, real const *extra_params) {

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

    RHS_cpu(sv, rDY, stim_curr, *dt, transmurality, sf_Iks, extra_params);
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
        RHS_cpu(sv, rDY, stim_curr, *dt, transmurality,sf_Iks, extra_params);
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

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real sf_Iks, real const *extra_params) {

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

    // State variables
    real v       = sv[0];
    real nai     = sv[1];
    real nass    = sv[2];
    real ki      = sv[3];
    real kss     = sv[4];
    real cai     = sv[5];
    real cass    = sv[6];
    real cansr   = sv[7];
    real cajsr   = sv[8];
    real m       = sv[9];
    real hp      = sv[10];
    real h       = sv[11];
    real j       = sv[12];
    real jp      = sv[13];
    real mL      = sv[14];
    real hL      = sv[15];
    real hLp     = sv[16];
    real a       = sv[17];
    real iF      = sv[18];
    real iS      = sv[19];
    real ap      = sv[20];
    real iFp     = sv[21];
    real iSp     = sv[22];
    real d       = sv[23];
    real ff      = sv[24];
    real fs      = sv[25];
    real fcaf    = sv[26];
    real fcas    = sv[27];
    real jca     = sv[28];
    real nca     = sv[29];
    real nca_i   = sv[30];
    real ffp     = sv[31];
    real fcafp   = sv[32];
    real xs1     = sv[33];
    real xs2     = sv[34];
    real Jrel_np = sv[35];
    real CaMKt   = sv[36];
    real ikr_c0  = sv[37];
    real ikr_c1  = sv[38];
    real ikr_c2  = sv[39];
    real ikr_o   = sv[40];
    real ikr_i   = sv[41];
    real Jrel_p  = sv[42];

    #include "ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.common.c"
}

void RHS_RL_cpu(real *a_, real *b_, const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real sf_Iks, real const *extra_params) {

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

    // State variables
    real v       = sv[0];
    real nai     = sv[1];
    real nass    = sv[2];
    real ki      = sv[3];
    real kss     = sv[4];
    real cai     = sv[5];
    real cass    = sv[6];
    real cansr   = sv[7];
    real cajsr   = sv[8];
    real m       = sv[9];
    real hp      = sv[10];
    real h       = sv[11];
    real j       = sv[12];
    real jp      = sv[13];
    real mL      = sv[14];
    real hL      = sv[15];
    real hLp     = sv[16];
    real a       = sv[17];
    real iF      = sv[18];
    real iS      = sv[19];
    real ap      = sv[20];
    real iFp     = sv[21];
    real iSp     = sv[22];
    real d       = sv[23];
    real ff      = sv[24];
    real fs      = sv[25];
    real fcaf    = sv[26];
    real fcas    = sv[27];
    real jca     = sv[28];
    real nca     = sv[29];
    real nca_i   = sv[30];
    real ffp     = sv[31];
    real fcafp   = sv[32];
    real xs1     = sv[33];
    real xs2     = sv[34];
    real Jrel_np = sv[35];
    real CaMKt   = sv[36];
    real ikr_c0  = sv[37];
    real ikr_c1  = sv[38];
    real ikr_c2  = sv[39];
    real ikr_o   = sv[40];
    real ikr_i   = sv[41];
    real Jrel_p  = sv[42];

    #include "ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments_RL.common.c"
}
