#include "trovato_2019.h"
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

    log_to_stdout_and_file("Using Trovato_2019 CPU model\n");

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
    for(uint32_t i = 0; i < num_cells; i++) 
    {
        real *sv = &solver->sv[i * NEQ];
        
        sv[0] = -86.6814002878592;
        sv[1] = 0.00505983330678751;
        sv[2] = 0.000101777993438818;
        sv[3] = 8.23183964616932;
        sv[4] = 8.23153516580562;
        sv[5] = 8.23154325237268;
        sv[6] = 143.767359809132;
        sv[7] = 143.767768218104;
        sv[8] = 143.767769906216;
        sv[9] = 4.36004404734282e-5;
        sv[10] = 0.000102004317781147;
        sv[11] = 1.26350902016858;
        sv[12] = 1.24811940209535;
        sv[13] = 1.26516959198518;
        sv[14] = 0.000108240945806962;
        sv[15] = 1.25045800437317e-69;
        sv[16] = 0.00632661703915808;
        sv[17] = 0.788611739889677;
        sv[18] = 0.788545979951331;
        sv[19] = 0.790474358603666;
        sv[20] = 0.579693514309867;
        sv[21] = 0.790947058236417;
        sv[22] = 0.000241925773627233;
        sv[23] = 0.463574582508218;
        sv[24] = 0.240216198686475;
        sv[25] = 0.000272851144435704;
        sv[26] = 0.649604795721571;
        sv[27] = 0.989965695822495;
        sv[28] = 6.97735089296892e-9;
        sv[29] = 0.999999968230738;
        sv[30] = 0.926692153319136;
        sv[31] = 0.99999996819573;
        sv[32] = 0.999999905741936;
        sv[33] = 0.999978907334662;
        sv[34] = 0.999999968365903;
        sv[35] = 0.999999968278239;
        sv[36] = 0.00547252500964926;
        sv[37] = 0.000304250912559619;
        sv[38] = 0.994214357917907;
        sv[39] = 0.000331691184084272;
        sv[40] = 0.568716473334161;
        sv[41] = 0.191165248085394;
        sv[42] = 0.000222677365291219;
        sv[43] = 0.233119011214908;
        sv[44] = 0.997084813729909;
        sv[45] = 0.466236137183558;
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
    for (uint32_t i = 0; i < num_cells_to_solve; i++) 
    {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) 
        {
            solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id);
        }
        else 
        {
            for (int j = 0; j < num_steps; ++j) 
            {
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

    // Forward Euler
    for (int i = 0; i < NEQ; i++)
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

    #include "trovato_common_inc.c"
}


