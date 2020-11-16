#include "Maleckar2009.h"
#include <stdlib.h>
real max_step;
real min_step;
real abstol;
real reltol;
bool adpt;
real *ode_dt, *ode_previous_dt, *ode_time_new;

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_to_stdout_and_file("Using Maleckar2009 CPU model\n");

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    max_step = solver->max_dt;
    min_step = solver->min_dt;
    abstol   = solver->abs_tol;
    reltol   = solver->rel_tol;
    adpt     = solver->adaptive;

    if(adpt) {
        ode_dt = (real*)malloc(num_cells*sizeof(real));

        OMP(parallel for)
        for(int i = 0; i < num_cells; i++) {
            ode_dt[i] = solver->min_dt;
        }

        ode_previous_dt = (real*)calloc(num_cells, sizeof(real));
        ode_time_new    = (real*)calloc(num_cells, sizeof(real));
        log_to_stdout_and_file("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_to_stdout_and_file("Using Euler model to solve the ODEs\n");
    }
    OMP(parallel for)
		for(uint32_t i = 0; i < num_cells; i++) {

			real *sv = &solver->sv[i * NEQ];

			sv[0] = -73.941851;
			sv[1] = 0.003325;
			sv[2] = 0.875262;
			sv[3] = 0.870692;
			sv[4] = 0.000014;
			sv[5] = 0.998578;
			sv[6] = 0.998561;
			sv[7] = 0.001098;
			sv[8] = 0.948202;
			sv[9] = 0.000371;
			sv[10] = 0.966869;
			sv[11] = 0.004661;
			sv[12] = 0.000054;
			sv[13] = 8.488527;
			sv[14] = 6.5e-5;
			sv[15] = 129.502075;
			sv[16] = 7.1e-5;
			sv[17] = 0.026604;
			sv[18] = 0.012843;
			sv[19] = 0.190077;
			sv[20] = 0.714719;
			sv[21] = 1.38222;
			sv[22] = 130.019282;
			sv[23] = 1.814418;
			sv[24] = 5.588239;
			sv[25] = 0.630471;
			sv[26] = 0.646226;
			sv[27] = 0.43071;
			sv[28] = 0.45453;
			sv[29] = 0.002665;

		}
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    #pragma omp parallel for private(sv_id)
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {

            solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id);
        }
        else {
            for (int j = 0; j < num_steps; ++j) {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);
            }

        }

    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id) {

    real rDY[NEQ];

    real _tolerances_[NEQ];
    real _aux_tol = 0.0;
    //initializes the variables
    real dt = ode_dt[sv_id];
    real time_new = ode_time_new[sv_id];
    real previous_dt = ode_previous_dt[sv_id];

    real edos_old_aux_[NEQ];
    real edos_new_euler_[NEQ];
    real *_k1__ = (real*) malloc(sizeof(real)*NEQ);
    real *_k2__ = (real*) malloc(sizeof(real)*NEQ);
    real *_k_aux__;

    const real _beta_safety_ = 0.8;

    const real __tiny_ = pow(abstol, 2.0f);

    if(time_new + dt > final_time) {
       dt = final_time - time_new;
    }

    RHS_cpu(sv, rDY, stim_curr, dt);
    time_new += dt;

    for(int i = 0; i < NEQ; i++){
        _k1__[i] = rDY[i];
    }

    int count = 0;

    int count_limit = (final_time - time_new)/min_step;

    int aux_count_limit = count_limit+2000000;

    if(aux_count_limit > 0) {
        count_limit = aux_count_limit;
    }

    while(1) {

        for(int i = 0; i < NEQ; i++) {
            //stores the old variables in a vector
            edos_old_aux_[i] = sv[i];
            //computes euler method
            edos_new_euler_[i] = _k1__[i] * dt + edos_old_aux_[i];
            //steps ahead to compute the rk2 method
            sv[i] = edos_new_euler_[i];
        }

        time_new += dt;
        RHS_cpu(sv, rDY, stim_curr, dt);
        time_new -= dt;//step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < NEQ; i++) {
            // stores the new evaluation
            _k2__[i] = rDY[i];
            _aux_tol = fabs(edos_new_euler_[i]) * reltol;
            _tolerances_[i] = (abstol > _aux_tol) ? abstol : _aux_tol;

            // finds the greatest error between  the steps
            auxError = fabs(((dt / 2.0) * (_k1__[i] - _k2__[i])) / _tolerances_[i]);

            greatestError = (auxError > greatestError) ? auxError : greatestError;
        }
        ///adapt the time step
        greatestError += __tiny_;
        previous_dt = dt;
        ///adapt the time step
        dt = _beta_safety_ * dt * sqrt(1.0f/greatestError);

        if (time_new + dt > final_time) {
            dt = final_time - time_new;
        }

        //it doesn't accept the solution
        if ( count < count_limit  && (greatestError >= 1.0f)) {
            //restore the old values to do it again
            for(int i = 0;  i < NEQ; i++) {
                sv[i] = edos_old_aux_[i];
            }
            count++;
            //throw the results away and compute again
        } else{//it accepts the solutions
            count = 0;

            if (dt < min_step) {
                dt = min_step;
            }

            else if (dt > max_step && max_step != 0) {
                dt = max_step;
            }

            if (time_new + dt > final_time) {
                dt = final_time - time_new;
            }

            _k_aux__ = _k2__;
            _k2__	= _k1__;
            _k1__	= _k_aux__;

            //it steps the method ahead, with euler solution
            for(int i = 0; i < NEQ; i++){
                sv[i] = edos_new_euler_[i];
            }

            if(time_new + previous_dt >= final_time) {
                if((fabs(final_time - time_new) < 1.0e-5)) {
                    break;
                } else if(time_new < final_time) {
                    dt = previous_dt = final_time - time_new;
                    time_new += previous_dt;
                    break;
                } else {
                    dt = previous_dt = min_step;
                    time_new += (final_time - time_new);
                    printf("Error: %lf\n", final_time - time_new);
                    break;
                }
            } else {
                time_new += previous_dt;
            }

        }
    }

	ode_dt[sv_id] = dt;
	ode_time_new[sv_id] = time_new;
	ode_previous_dt[sv_id] = previous_dt;
	
    free(_k1__);
    free(_k2__);
}

void RHS_cpu(const real *sv, real *rDY, real stim_current, real dt) {

	real var_membrane__V = sv[0]; // Units: millivolt; Initial value: -73.941851
	real var_sodium_current_m_gate__m = sv[1]; // Units: dimensionless; Initial value: 0.003325
	real var_sodium_current_h1_gate__h1 = sv[2]; // Units: dimensionless; Initial value: 0.875262
	real var_sodium_current_h2_gate__h2 = sv[3]; // Units: dimensionless; Initial value: 0.870692
	real var_L_type_Ca_channel_d_L_gate__d_L = sv[4]; // Units: dimensionless; Initial value: 0.000014
	real var_L_type_Ca_channel_f_L1_gate__f_L1 = sv[5]; // Units: dimensionless; Initial value: 0.998578
	real var_L_type_Ca_channel_f_L2_gate__f_L2 = sv[6]; // Units: dimensionless; Initial value: 0.998561
	real var_Ca_independent_transient_outward_K_current_r_gate__r = sv[7]; // Units: dimensionless; Initial value: 0.001098
	real var_Ca_independent_transient_outward_K_current_s_gate__s = sv[8]; // Units: dimensionless; Initial value: 0.948202
	real var_ultra_rapid_K_current_aur_gate__a_ur = sv[9]; // Units: dimensionless; Initial value: 0.000371
	real var_ultra_rapid_K_current_iur_gate__i_ur = sv[10]; // Units: dimensionless; Initial value: 0.966869
	real var_delayed_rectifier_K_currents_n_gate__n = sv[11]; // Units: dimensionless; Initial value: 0.004661
	real var_delayed_rectifier_K_currents_pa_gate__pa = sv[12]; // Units: dimensionless; Initial value: 0.000054
	real var_intracellular_ion_concentrations__Na_i = sv[13]; // Units: millimolar; Initial value: 8.488527
	real var_intracellular_ion_concentrations__Ca_i = sv[14]; // Units: millimolar; Initial value: 6.5e-5
	real var_intracellular_ion_concentrations__K_i = sv[15]; // Units: millimolar; Initial value: 129.502075
	real var_intracellular_ion_concentrations__Ca_d = sv[16]; // Units: millimolar; Initial value: 7.1e-5
	real var_intracellular_Ca_buffering__O_C = sv[17]; // Units: dimensionless; Initial value: 0.026604
	real var_intracellular_Ca_buffering__O_TC = sv[18]; // Units: dimensionless; Initial value: 0.012843
	real var_intracellular_Ca_buffering__O_TMgC = sv[19]; // Units: dimensionless; Initial value: 0.190077
	real var_intracellular_Ca_buffering__O_TMgMg = sv[20]; // Units: dimensionless; Initial value: 0.714719
	real var_cleft_space_ion_concentrations__Na_c = sv[22]; // Units: millimolar; Initial value: 130.019282
	real var_cleft_space_ion_concentrations__Ca_c = sv[23]; // Units: millimolar; Initial value: 1.814418
	real var_cleft_space_ion_concentrations__K_c = sv[24]; // Units: millimolar; Initial value: 5.588239
	real var_Ca_handling_by_the_SR__Ca_rel = sv[25]; // Units: millimolar; Initial value: 0.630471
	real var_Ca_handling_by_the_SR__Ca_up = sv[26]; // Units: millimolar; Initial value: 0.646226
	real var_Ca_handling_by_the_SR__O_Calse = sv[27]; // Units: dimensionless; Initial value: 0.43071
	real var_Ca_handling_by_the_SR__F1 = sv[28]; // Units: dimensionless; Initial value: 0.45453
	real var_Ca_handling_by_the_SR__F2 = sv[29]; // Units: dimensionless; Initial value: 0.002665	

	#include "Maleckar2009_common.inc.c"

}
