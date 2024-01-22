#include <stdio.h>
#include "CRN2017_TASK1_IKACh.h"

real max_step;
real min_step;
real abstol;
real reltol;
bool adpt;
real *ode_dt, *ode_previous_dt, *ode_time_new;

GET_CELL_MODEL_DATA(init_cell_model_data) 
{

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) 
{

    log_info("Using CRN2017 Regions CPU model\n");

    uint32_t num_cells = solver->original_num_cells;
    
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    max_step = solver->max_dt;
    min_step = solver->min_dt;
    abstol   = solver->abs_tol;
    reltol   = solver->rel_tol;
    adpt     = solver->adaptive;

    if(adpt) 
    {
        ode_dt = (real*)malloc(num_cells*sizeof(real));

        OMP(parallel for)
        for(int i = 0; i < num_cells; i++) 
        {
            ode_dt[i] = solver->min_dt;
        }

        ode_previous_dt = (real*)calloc(num_cells, sizeof(real));
        ode_time_new    = (real*)calloc(num_cells, sizeof(real));
        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } 
    else 
    {
        log_info("Using Euler model to solve the ODEs\n");
    }

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) 
    {

        real *sv = &solver->sv[i * NEQ];
        if(solver->ode_extra_data == NULL) {
            // Initial conditions CRN 2017 version
            sv[0]       = -82.7618;     
            sv[1]       =   0.0022;
            sv[2]       =   0.9752;
            sv[3]       =   0.9844;
            sv[4]       =   0.0279;
            sv[5]       =   0.9994;
            sv[6]       =   1.1219e-4;
            sv[7]       =   0.9070;
            sv[8]       =   0.7025;
            sv[9]       =   0.0042;
            sv[10]      =   0.7816;
            sv[11]      =   0.9488;
            sv[12]      =   0.0044;
            sv[13]      =   0.0184;
            sv[14]      =   0.0;
            // Ca handling
            sv[15]      =   0.7719;      // mM 
            sv[16]      =   1.1355;      // mM
            sv[17]      =   0.0000;
            sv[18]      =   1.0000;
            sv[19]      =   0.9993;
            // Concentrations
            sv[20]      =  11.9137;
            sv[21]      = 140.0000;
            sv[22]      = 138.0694; 
            sv[23]      =   5.4000;     // mM
            sv[24]      =   1.4810e-04; // mM
            sv[25]      =   1.8000;     // mM 
        }
        else {
            real *extra_parameters = (real *)solver->ode_extra_data;
            sv[ 0] = extra_parameters[14];
            sv[ 1] = extra_parameters[15];
            sv[ 2] = extra_parameters[16];
            sv[ 3] = extra_parameters[17];
            sv[ 4] = extra_parameters[18];
            sv[ 5] = extra_parameters[19];
            sv[ 6] = extra_parameters[20];
            sv[ 7] = extra_parameters[21];
            sv[ 8] = extra_parameters[22];
            sv[ 9] = extra_parameters[23];
            sv[10] = extra_parameters[24];
            sv[11] = extra_parameters[25];
            sv[12] = extra_parameters[26];
            sv[13] = extra_parameters[27];
            sv[14] = extra_parameters[28];
            sv[15] = extra_parameters[29];
            sv[16] = extra_parameters[30];
            sv[17] = extra_parameters[31];
            sv[18] = extra_parameters[32];
            sv[19] = extra_parameters[33];
            sv[20] = extra_parameters[34];
            sv[21] = extra_parameters[35];
            sv[22] = extra_parameters[36];
            sv[23] = extra_parameters[37];
            sv[24] = extra_parameters[38];
            sv[25] = extra_parameters[39];
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

    // Get the extra data array
    int offset = 40;
    real *extra_parameters = NULL;
    if(ode_solver->ode_extra_data) 
    {
        extra_parameters = ((real*)ode_solver->ode_extra_data);
    }
    else 
    {
        log_error_and_exit("You need to specify a mask function when using this mixed model!\n");
    }

    #pragma omp parallel for private(sv_id)
    for (int i = 0; i < num_cells_to_solve; i++) 
    {
        int mapping = extra_parameters[i+offset];

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) 
        {
            solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver->ode_extra_data, mapping);
        }
        else 
        {
            for (int j = 0; j < num_steps; ++j) 
            {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], ode_solver->ode_extra_data, mapping);
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current, real *extra_parameters, int mapping)  
{

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, extra_parameters, mapping);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, real *extra_parameters, int mapping) 
{

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

    RHS_cpu(sv, rDY, stim_curr, dt, extra_parameters, mapping);
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
        RHS_cpu(sv, rDY, stim_curr, dt, extra_parameters, mapping);
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
            _k2__   = _k1__;
            _k1__   = _k_aux__;

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

void RHS_cpu(const real *sv, real *rDY, real stim_current, real dt, real *extra_parameters, int mapping) {

    // State variables
    real Vm         = sv[ 0];
    real INa_va     = sv[ 1];
    real INa_vi_1   = sv[ 2];
    real INa_vi_2   = sv[ 3];
    real Ito_va     = sv[ 4];
    real Ito_vi     = sv[ 5];
    real ICaL_va    = sv[ 6];
    real ICaL_vi    = sv[ 7];
    real ICaL_ci    = sv[ 8];
    real IKur_va    = sv[ 9];
    real IKur_viS   = sv[10];
    real IKur_viF   = sv[11];
    real IKr_va     = sv[12];
    real IKs_va     = sv[13];
    real IK2P_va    = sv[14];
    real CajSR      = sv[15];
    real CanSR      = sv[16];
    real RyRo       = sv[17];
    real RyRr       = sv[18];
    real RyRi       = sv[19];
    real Nai        = sv[20];
    real Nao        = sv[21];
    real Ki         = sv[22];
    real Ko         = sv[23];
    real Cai        = sv[24];
    real Cao        = sv[25];

    #include "CRN2017_TASK1_IKACh_common.inc.c"

}

