#include <assert.h>
#include <stdlib.h>
#include "fhn_mitchell.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V_1;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ_1;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    static bool first_call = true;

    if(first_call) {
        log_to_stdout_and_file("Using mixed version of modified FHN 1961 + Mitchell-Shaeffer 2003 CPU model\n");
        first_call = false;
    }
    
    // Get the mapping array
    uint32_t *mapping = NULL;
    if(ode_extra_config) {
        mapping = (uint32_t*)(solver->ode_extra_data);
    }
    else {
        log_to_stderr_and_file_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    uint32_t num_volumes = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ_1*num_volumes*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_volumes; i++) {

        // FHN
        if (mapping[i] == 0) {
            real *sv = &solver->sv[i * NEQ_1];
            sv[0] = 0.0;    // Vm millivolt
            sv[1] = 0.0;    // h dimensionless
        }
        // Mitchell-Shaeffer
        else if (mapping[i] == 1) {
            real *sv = &solver->sv[i * NEQ_2];
            sv[0] = 0.00000820413566106744;     // Vm millivolt
            sv[1] = 0.8789655121804799f;        // h dimensionless
        }
        else {
            log_to_stderr_and_file_and_exit("[ERR] Invalid mapping value '%u'! It should be 0(FHN) or 1(Mitchell-Shaeffer)!\n",mapping[i]);
        }
    }

}

SOLVE_MODEL_ODES (solve_model_odes_cpu) {

    // Get the mapping array
    uint32_t *mapping = NULL;
    if(ode_extra_config) {
        mapping = (uint32_t*)(ode_solver->ode_extra_data);
    }
    else {
        log_to_stderr_and_file_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (uint32_t j = 0; j < num_steps; ++j) {

            // FHN
            if (mapping[i] == 0) {
                solve_model_ode_cpu_fhn(dt, sv + (sv_id * NEQ_1), stim_currents[i]);
            }
            // Mitchell-Shaeffer
            else if (mapping[i] == 1) {
                solve_model_ode_cpu_mitchell(dt, sv + (sv_id * NEQ_2), stim_currents[i]);
            }
            else {
                log_to_stderr_and_file_and_exit("[ERR] Invalid mapping value '%u'! It should be 0(FHN) or 1(Mitchell-Shaeffer)!\n",mapping[i]);
            }
            
        }
    }

}

void solve_model_ode_cpu_fhn (real dt, real *sv, real stim_current)  {

    real rY[NEQ_1], rDY[NEQ_1];

    for(int i = 0; i < NEQ_1; i++)
        rY[i] = sv[i];

    RHS_cpu_fhn(rY, rDY, stim_current);

    for(int i = 0; i < NEQ_1; i++)
        sv[i] = dt*rDY[i] + rY[i];
    
}

void RHS_cpu_fhn (const real *sv, real *rDY_, real stim_current) {

    //State variables
    const real u = sv[0];
    const real v = sv[1];

    // Constants
    const real a = 0.2f;
    const real b = 0.5f;
    const real k = 36.0;
    const real epsilon  =  0.000150;

    // Rates
    rDY_[0] = k*(u*(1.0f - u)*(u - a) - u*v) + stim_current;
    rDY_[1] = k*epsilon*(b*u - v);

}

void solve_model_ode_cpu_mitchell (real dt, real *sv, real stim_current) {
    real rY[NEQ_2], rDY[NEQ_2];

    for(int i = 0; i < NEQ_2; i++)
        rY[i] = sv[i];

    RHS_cpu_mitchell(rY, rDY, stim_current);

    for(int i = 0; i < NEQ_2; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu_mitchell(const real *sv, real *rDY_, real stim_current) {
    //State variables
    const real V = sv[0];
    const real h = sv[1];

    // Constants
    const real tau_in = 0.3;
    const real tau_out = 6.0;
    const real V_gate = 0.13;
    const real tau_open = 120.0;
    const real tau_close = 150.0;

    // Algebraics
    real J_stim = stim_current;
    real J_in = ( h*( pow(V, 2.00000)*(1.00000 - V)))/tau_in;
    real J_out = - (V/tau_out);

    // Rates
    rDY_[0] = J_out + J_in + J_stim;
    rDY_[1] = (V < V_gate ? (1.00000 - h)/tau_open : - h/tau_close);
}