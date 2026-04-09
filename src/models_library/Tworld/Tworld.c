#include "Tworld.h"
#include <stdlib.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
    cell_model->initial_v = INITIAL_V;
    if(get_neq)
    cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_msg("Using Tworld CPU model\n");

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));
    bool adpt = solver->adaptive;

    if (adpt) {
        solver->ode_dt = (real*)malloc(num_cells*sizeof(real));

        OMP(parallel for)
        for(int i = 0; i < num_cells; i++) {
            solver->ode_dt[i] = solver->min_dt;
        }

        solver->ode_previous_dt = (real*)calloc(num_cells, sizeof(real));
        solver->ode_time_new    = (real*)calloc(num_cells, sizeof(real));
        log_info("Using Adaptive timestep model to solve the ODEs\n");
    }
    else {
        log_info("Using Fixed timestep to solve the ODEs\n");
    }

    real *extra_data = NULL;
    if(solver->ode_extra_data) {
        extra_data = (real*)solver->ode_extra_data;
    } else {
        log_error_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {
        real *sv = &solver->sv[i * NEQ];
		int layer = (int) extra_data[i];

        #include "TWorld_human_V9_SS.common.c"
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
	
    real *extra_data = NULL;
    if(ode_solver->ode_extra_data) {
        extra_data = (real*)ode_solver->ode_extra_data;
    } else {
        log_error_and_exit("You need to specify a mask function when using this mixed model!\n");
    }


    #pragma omp parallel for private(sv_id)
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {
		
        int layer = (int) 			extra_data[i];
       	real transmural = 			extra_data[i + (1*num_cells_to_solve)];
		real apicobasal = 			extra_data[i + (2*num_cells_to_solve)];
		real ischemia = 			extra_data[i + (3*num_cells_to_solve)];
		
		int modelVF	 = (int)			extra_data[7*num_cells_to_solve + 0];
		int NZ_severity = (int) 		extra_data[7*num_cells_to_solve + 1];
		int ICZ_severity = (int) 		extra_data[7*num_cells_to_solve + 2];		

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver, 
								layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);
        }
        else {
            for (int j = 0; j < num_steps; ++j) {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], 
											layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current, int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity) {

    const real TOLERANCE = 1e-08;
    real rY[NEQ], rDY[NEQ];
    real a[NEQ], b[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_RL_cpu(a, b, rY, rDY, stim_current, dt, 
					layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);

    // Solve variables based on its type:
    //  Non-linear = Euler
    //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
    for (int i = 0; i <= 23; i++)
        SOLVE_EQUATION_EULER_CPU(i);
    for (int i = 24; i <= 39; i++)
        sv[i] = rDY[i];
    for (int i = 40; i <= 44; i++)
        SOLVE_EQUATION_EULER_CPU(i);
    for (int i = 45; i <= 50; i++)
        sv[i] = rDY[i];
    for (int i = 51; i <= 52; i++)
        SOLVE_EQUATION_EULER_CPU(i);
    for (int i = 53; i <= 54; i++)
        sv[i] = rDY[i];
    for (int i = 55; i <= 64; i++)
        SOLVE_EQUATION_EULER_CPU(i);
    for (int i = 65; i <= 67; i++)
        sv[i] = rDY[i];
    for (int i = 68; i <= 74; i++)
        SOLVE_EQUATION_EULER_CPU(i);
    for (int i = 75; i <= 77; i++)
        sv[i] = rDY[i];
    SOLVE_EQUATION_EULER_CPU(78);
    for (int i = 79; i <= 92; i++)
        sv[i] = rDY[i];

}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver, 
									int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity) {

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
   
    RHS_cpu(sv, rDY, stim_curr, *dt, 
						layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);
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
        RHS_cpu(sv, rDY, stim_curr, *dt, 
						layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);
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

void RHS_cpu (const real *sv, real *rDY, real stim_current, real dt, 
				int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity){

    // State variables
    real_cpu v = sv[0];
    real_cpu ions_na_junc = sv[1];
    real_cpu ions_na_sl = sv[2];
    real_cpu ions_na_i = sv[3];
    real_cpu ions_k_i = sv[4];
    real_cpu ions_ca_junc = sv[5];
    real_cpu ions_ca_sl = sv[6];
    real_cpu ions_ca_i = sv[7];
    real_cpu ions_cl_i = sv[8];
    real_cpu ions_ca_SR = sv[9];
    real_cpu buffers_NaBj = sv[10];
    real_cpu buffers_NaBsl = sv[11];
    real_cpu buffers_TnClow = sv[12];
    real_cpu buffers_TnCHc = sv[13];
    real_cpu buffers_TnCHm = sv[14];
    real_cpu buffers_CaM = sv[15];
    real_cpu buffers_Myosin_ca = sv[16];
    real_cpu buffers_Myosin_mg = sv[17];
    real_cpu buffers_SRB = sv[18];
    real_cpu buffers_SLLj = sv[19];
    real_cpu buffers_SLLsl = sv[20];
    real_cpu buffers_SLHj = sv[21];
    real_cpu buffers_SLHsl = sv[22];
    real_cpu buffers_Csqn = sv[23];
    real_cpu m = sv[24];
    real_cpu h = sv[25];
    real_cpu j = sv[26];
    real_cpu mL = sv[27];
    real_cpu hL = sv[28];
    real_cpu hLp = sv[29];
    real_cpu ito_xtos = sv[30];
    real_cpu ito_ytos = sv[31];
    real_cpu ito_xtof = sv[32];
    real_cpu ito_ytof = sv[33];
    real_cpu ito_xtos_p = sv[34];
    real_cpu ito_ytos_p = sv[35];
    real_cpu ito_xtof_p = sv[36];
    real_cpu ito_ytof_p = sv[37];
    real_cpu iks_xs_junc = sv[38];
    real_cpu iks_xs_sl = sv[39];
    real_cpu c0 = sv[40];
    real_cpu c1 = sv[41];
    real_cpu c2 = sv[42];
    real_cpu o = sv[43];
    real_cpu i = sv[44];
    real_cpu d = sv[45];
    real_cpu ff = sv[46];
    real_cpu fs = sv[47];
    real_cpu fcaf = sv[48];
    real_cpu fcas = sv[49];
    real_cpu jca = sv[50];
    real_cpu nca = sv[51];
    real_cpu nca_i = sv[52];
    real_cpu ffp = sv[53];
    real_cpu fcafp = sv[54];
    real_cpu ical_pureCDI_junc = sv[55];
    real_cpu ical_pureCDI_sl = sv[56];
    real_cpu ryr_R = sv[57];
    real_cpu ryr_O = sv[58];
    real_cpu ryr_I = sv[59];
    real_cpu ryr_CaRI = sv[60];
    real_cpu ryr_R_p = sv[61];
    real_cpu ryr_O_p = sv[62];
    real_cpu ryr_I_p = sv[63];
    real_cpu ryr_CaRI_p = sv[64];
    real_cpu jrel_icaldep_act = sv[65];
    real_cpu jrel_icaldep_f1 = sv[66];
    real_cpu jrel_icaldep_f2 = sv[67];
    real_cpu contraction_XS = sv[68];
    real_cpu contraction_XW = sv[69];
    real_cpu contraction_Ca_TRPN = sv[70];
    real_cpu contraction_TmBlocked = sv[71];
    real_cpu contraction_ZETAS = sv[72];
    real_cpu contraction_ZETAW = sv[73];
    real_cpu camk_trap = sv[74];
    real_cpu camk_f_ICaL = sv[75];
    real_cpu camk_f_RyR = sv[76];
    real_cpu camk_f_PLB = sv[77];
    real_cpu casig_serca_trap = sv[78];
    real_cpu hp = sv[79];
    real_cpu jp = sv[80];
    real_cpu m_P = sv[81];
    real_cpu h_P = sv[82];
    real_cpu j_P = sv[83];
    real_cpu hp_P = sv[84];
    real_cpu jp_P = sv[85];
    real_cpu d_P = sv[86];
    real_cpu ff_P = sv[87];
    real_cpu fs_P = sv[88];
    real_cpu fcaf_P = sv[89];
    real_cpu fcas_P = sv[90];
    real_cpu fBPf = sv[91];
    real_cpu fcaBPf = sv[92];

    #include "Tworld_common.inc.c"
}

void RHS_RL_cpu(real *a_, real *b_, const real *sv, real *rDY, real stim_current, real dt, int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity) 
{
    // State variables
    real_cpu v = sv[0];
    real_cpu ions_na_junc = sv[1];
    real_cpu ions_na_sl = sv[2];
    real_cpu ions_na_i = sv[3];
    real_cpu ions_k_i = sv[4];
    real_cpu ions_ca_junc = sv[5];
    real_cpu ions_ca_sl = sv[6];
    real_cpu ions_ca_i = sv[7];
    real_cpu ions_cl_i = sv[8];
    real_cpu ions_ca_SR = sv[9];
    real_cpu buffers_NaBj = sv[10];
    real_cpu buffers_NaBsl = sv[11];
    real_cpu buffers_TnClow = sv[12];
    real_cpu buffers_TnCHc = sv[13];
    real_cpu buffers_TnCHm = sv[14];
    real_cpu buffers_CaM = sv[15];
    real_cpu buffers_Myosin_ca = sv[16];
    real_cpu buffers_Myosin_mg = sv[17];
    real_cpu buffers_SRB = sv[18];
    real_cpu buffers_SLLj = sv[19];
    real_cpu buffers_SLLsl = sv[20];
    real_cpu buffers_SLHj = sv[21];
    real_cpu buffers_SLHsl = sv[22];
    real_cpu buffers_Csqn = sv[23];
    real_cpu m = sv[24];
    real_cpu h = sv[25];
    real_cpu j = sv[26];
    real_cpu mL = sv[27];
    real_cpu hL = sv[28];
    real_cpu hLp = sv[29];
    real_cpu ito_xtos = sv[30];
    real_cpu ito_ytos = sv[31];
    real_cpu ito_xtof = sv[32];
    real_cpu ito_ytof = sv[33];
    real_cpu ito_xtos_p = sv[34];
    real_cpu ito_ytos_p = sv[35];
    real_cpu ito_xtof_p = sv[36];
    real_cpu ito_ytof_p = sv[37];
    real_cpu iks_xs_junc = sv[38];
    real_cpu iks_xs_sl = sv[39];
    real_cpu c0 = sv[40];
    real_cpu c1 = sv[41];
    real_cpu c2 = sv[42];
    real_cpu o = sv[43];
    real_cpu i = sv[44];
    real_cpu d = sv[45];
    real_cpu ff = sv[46];
    real_cpu fs = sv[47];
    real_cpu fcaf = sv[48];
    real_cpu fcas = sv[49];
    real_cpu jca = sv[50];
    real_cpu nca = sv[51];
    real_cpu nca_i = sv[52];
    real_cpu ffp = sv[53];
    real_cpu fcafp = sv[54];
    real_cpu ical_pureCDI_junc = sv[55];
    real_cpu ical_pureCDI_sl = sv[56];
    real_cpu ryr_R = sv[57];
    real_cpu ryr_O = sv[58];
    real_cpu ryr_I = sv[59];
    real_cpu ryr_CaRI = sv[60];
    real_cpu ryr_R_p = sv[61];
    real_cpu ryr_O_p = sv[62];
    real_cpu ryr_I_p = sv[63];
    real_cpu ryr_CaRI_p = sv[64];
    real_cpu jrel_icaldep_act = sv[65];
    real_cpu jrel_icaldep_f1 = sv[66];
    real_cpu jrel_icaldep_f2 = sv[67];
    real_cpu contraction_XS = sv[68];
    real_cpu contraction_XW = sv[69];
    real_cpu contraction_Ca_TRPN = sv[70];
    real_cpu contraction_TmBlocked = sv[71];
    real_cpu contraction_ZETAS = sv[72];
    real_cpu contraction_ZETAW = sv[73];
    real_cpu camk_trap = sv[74];
    real_cpu camk_f_ICaL = sv[75];
    real_cpu camk_f_RyR = sv[76];
    real_cpu camk_f_PLB = sv[77];
    real_cpu casig_serca_trap = sv[78];
    real_cpu hp = sv[79];
    real_cpu jp = sv[80];
    real_cpu m_P = sv[81];
    real_cpu h_P = sv[82];
    real_cpu j_P = sv[83];
    real_cpu hp_P = sv[84];
    real_cpu jp_P = sv[85];
    real_cpu d_P = sv[86];
    real_cpu ff_P = sv[87];
    real_cpu fs_P = sv[88];
    real_cpu fcaf_P = sv[89];
    real_cpu fcas_P = sv[90];
    real_cpu fBPf = sv[91];
    real_cpu fcaBPf = sv[92];

    #include "Tworld_RL_common.inc.c"
}
