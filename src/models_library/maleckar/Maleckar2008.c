#include "Maleckar2008.h"
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

    log_info("Using Maleckar2008 CPU model\n");

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
        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_info("Using Euler model to solve the ODEs\n");
    }
    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {

        real *sv = &solver->sv[i * NEQ];

        sv[0]  = -87.169816169406;
        sv[1]  = 0.001075453357;
        sv[2]  = 0.990691306716;
        sv[3]  = 0.993888937283;
        sv[4]  = 0.000018211252;
        sv[5]  = 0.979322592773;
        sv[6]  = 0.001208153482;
        sv[7]  = 0.000033616596;
        sv[8]  = 0.004173008466;
        sv[9]  = 0.015242594688;
        sv[10] = 0.007074239331;
        sv[11] = 0.048267587131;
        sv[12] = 0.105468807033;
        sv[13] = 0.00364776906;
        sv[14] = 0.174403618112;
        sv[15] = 0.003643592594;
        sv[16] = 0.993331326442;
        sv[17] = 97.505463697266;
        sv[18] = 0.006679257264;
        sv[19] = 11.441712311614;
        sv[20] = 1.716573130685;
        sv[21] = 0.226941113355;
        sv[22] = 0.256752008084;
        sv[23] = 104.450004990523;
        sv[24] = 22.171689894953;
        sv[25] = 19.864701949854;
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



float __agos_factorial(int f){
	if(f>=0 && f<2)
		return 1.0;
	else if(f < 0)
		return 0.0/0.0;
	for(int i=f-1; i>=2; i--)
		f *= i;
	return (float)f;
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

    real var_cell__V = sv[0];
    real var_INa__xm = sv[1];
    real var_INa__xh = sv[2];
    real var_INa__xj = sv[3];
    real var_ICaL__c1 = sv[4];
    real var_ICaL__c2 = sv[5];
    real var_ICaL__xi1ca = sv[6];
    real var_ICaL__xi1ba = sv[7];
    real var_ICaL__xi2ca = sv[8];
    real var_ICaL__xi2ba = sv[9];
    real var_IKr__xr = sv[10];
    real var_IKs__xs1 = sv[11];
    real var_IKs__xs2 = sv[12];
    real var_Ito__xtos = sv[13];
    real var_Ito__ytos = sv[14];
    real var_Ito__xtof = sv[15];
    real var_Ito__ytof = sv[16];
    real var_Irel__Ca_JSR = sv[17];
    real var_Irel__xir = sv[18];
    real var_Na__Na_i = sv[19];
    real var_Ca__Ca_dyad = sv[20];
    real var_Ca__Ca_submem = sv[21];
    real var_Ca__Ca_i = sv[22];
    real var_Ca__Ca_NSR = sv[23];
    real var_Ca__tropi = sv[24];
    real var_Ca__trops = sv[25];

    #include "Maleckar2008_common.inc.c"

}

