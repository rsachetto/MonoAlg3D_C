// TODO: Describe all the variable names

#include "stewart_aslanidi_noble_2009.h"
#include <stdio.h>

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

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) 
{

    log_to_stdout_and_file("Using Stewart-Aslanidi-Noble 2009 CPU model\n");

    log_to_stdout_and_file("Using Stewart-Aslanidi-Noble 2009 CPU model\n");

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

        // Initial conditions from the original paper
        /*
        sv[0] = -74.7890522727;
        sv[1] = 136.9896086978;
        sv[2] = 8.5447311020;
        sv[3] = 0.0001720623;
        sv[4] = 0.0184308075;
        sv[5] = 0.4663168269;
        sv[6] = 0.3657472179;
        sv[7] = 0.0486609588;
        sv[8] = 0.0145766758;
        sv[9] = 0.2979720207;
        sv[10] = 0.0692509548;
        sv[11] = 0.0006146554;
        sv[12] = 0.0001356656;
        sv[13] = 0.5943228461;
        sv[14] = 0.8265709174;
        sv[15] = 0.9767040566;
        sv[16] = 0.9717098312;
        sv[17] = 0.0006830833;
        sv[18] = 3.2830723338;
        sv[19] = 0.8199969443;
        */

        // Steady-State for BCL=1000ms
        sv[0] = -69.1370441635924;
        sv[1] = 136.781894160227;
        sv[2] = 8.80420286531673;
        sv[3] = 0.000101878186157052;
        sv[4] = 0.0457562667986602;
        sv[5] = 0.00550281999719088;
        sv[6] = 0.313213286437995;
        sv[7] = 0.00953708522974789;
        sv[8] = 0.0417391656294997;
        sv[9] = 0.190678733735145;
        sv[10] = 0.238219836154029;
        sv[11] = 0.000446818714055411;
        sv[12] = 0.000287906256206415;
        sv[13] = 0.989328560287987;
        sv[14] = 0.995474890442185;
        sv[15] = 0.999955429598213;
        sv[16] = 0.96386101799501;
        sv[17] = 0.00103618091196912;
        sv[18] = 3.10836886659417;
        sv[19] = 0.991580051907845;

    }

}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    uint32_t sv_id;

    OMP(parallel for private(sv_id))
    for (uint32_t i = 0; i < num_cells_to_solve; i++) {

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

void solve_model_ode_cpu(real dt, real *sv, real stim_current) {

    real rY[NEQ], rDY[NEQ];

    // Save old value of the state vector
    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    // Solve Right-hand-side of the ODE's
    RHS_cpu(rY, rDY, stim_current, dt);

    // Solve model using Forward Euler
    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id) {

    const real _beta_safety_ = 0.8;
    int numEDO = NEQ;

    real rDY[numEDO];

    real _tolerances_[numEDO];
    real _aux_tol = 0.0;
    //initializes the variables
    ode_previous_dt[sv_id] = ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real*) malloc(sizeof(real)*numEDO);
    real *_k2__ = (real*) malloc(sizeof(real)*numEDO);
    real *_k_aux__;

    real *dt = &ode_dt[sv_id];
    real *time_new = &ode_time_new[sv_id];
    real *previous_dt = &ode_previous_dt[sv_id];

    if(*time_new + *dt > final_time) {
       *dt = final_time - *time_new;
    }

    RHS_cpu(sv, rDY, stim_curr, *dt);
    *time_new += *dt;

    for(int i = 0; i < numEDO; i++){
        _k1__[i] = rDY[i];
    }

    const double __tiny_ = pow(abstol, 2.0);

    int count = 0;

    int count_limit = (final_time - *time_new)/min_step;

    int aux_count_limit = count_limit+2000000;

    if(aux_count_limit > 0) {
        count_limit = aux_count_limit;
    }

    while(1) {

        for(int i = 0; i < numEDO; i++) {
            //stores the old variables in a vector
            edos_old_aux_[i] = sv[i];
            //computes euler method
            edos_new_euler_[i] = _k1__[i] * *dt + edos_old_aux_[i];
            //steps ahead to compute the rk2 method
            sv[i] = edos_new_euler_[i];
        }

        *time_new += *dt;
        RHS_cpu(sv, rDY, stim_curr, *dt);
        *time_new -= *dt;//step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < numEDO; i++) {
            //stores the new evaluation
            _k2__[i] = rDY[i];
            _aux_tol = fabs(edos_new_euler_[i])*reltol;
            _tolerances_[i] = (abstol > _aux_tol )?abstol:_aux_tol;
            //finds the greatest error between  the steps
            auxError = fabs(( (*dt/2.0)*(_k1__[i] - _k2__[i])) / _tolerances_[i]);

            greatestError = (auxError > greatestError) ? auxError : greatestError;
        }
        ///adapt the time step
        greatestError += __tiny_;
        *previous_dt = *dt;
        ///adapt the time step
        *dt = _beta_safety_ * (*dt) * sqrt(1.0f/greatestError);

        if (*time_new + *dt > final_time) {
            *dt = final_time - *time_new;
        }

        //it doesn't accept the solution
        if ( count < count_limit  && (greatestError >= 1.0f)) {
            //restore the old values to do it again
            for(int i = 0;  i < numEDO; i++) {
                sv[i] = edos_old_aux_[i];
            }

            count++;
            //throw the results away and compute again
        } else{//it accepts the solutions


            if(greatestError >=1.0) {
                printf("Accepting solution with error > %lf \n", greatestError);
            }

            //printf("%e %e\n", _ode->time_new, edos_new_euler_[0]);
            if (*dt < min_step) {
                *dt = min_step;
            }

            else if (*dt > max_step && max_step != 0) {
                *dt = max_step;
            }

            if (*time_new + *dt > final_time) {
                *dt = final_time - *time_new;
            }

            _k_aux__ = _k2__;
            _k2__	= _k1__;
            _k1__	= _k_aux__;

            //it steps the method ahead, with euler solution
            for(int i = 0; i < numEDO; i++){
                sv[i] = edos_new_euler_[i];
            }

            if(*time_new + *previous_dt >= final_time){
                if((fabs(final_time - *time_new) < 1.0e-5) ){
                    break;
                }else if(*time_new < final_time){
                    *dt = *previous_dt = final_time - *time_new;
                    *time_new += *previous_dt;
                    break;

                }else{
                    printf("Error: time_new %.20lf final_time %.20lf diff %e \n", *time_new , final_time, fabs(final_time - *time_new) );
                    break;
                }
            }else{
                *time_new += *previous_dt;
            }

        }
    }

    free(_k1__);
    free(_k2__);
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    real STATES[NEQ];
    for (uint32_t i = 0; i < NEQ; i++)
        STATES[i] = sv[i]; 

    #include "stewart_aslanidi_noble_2009_common.inc"
}

// The automatic pacing from the Purkinje cells can be interrupted by blocking the INa current by 100% 
// (ALGEBRAIC[54] = INa)
