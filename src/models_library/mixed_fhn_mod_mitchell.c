#include <stdio.h>
#include "mixed_fhn_mod_mitchell.h"

// TODO: Maybe change this function
//       Set number_of_ode_equations to the maximum 'NEQ' ?
GET_CELL_MODEL_DATA(init_cell_model_data) 
{

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V_1;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ_1;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) 
{

    static bool first_call = true;

    if(first_call) 
    {
        print_to_stdout_and_file("Using mixed version of modified FHN 1961 + Mitchell-Shaeffer 2003 CPU model\n");
        first_call = false;
    }
    
    // Get the mapping array
    uint32_t *mapping = NULL;
    if(extra_data) 
    {
        mapping = (uint32_t*)extra_data;
    }
    else 
    {
        print_to_stderr_and_file_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    // Based on the mapping initialize the initial conditions from the correct celular model
    if (mapping[sv_id] == 0)
    {
        sv[0] = 0.000000f; //Vm millivolt
        sv[1] = 0.000000f; //v dimensionless
    }
    else
    {
        sv[0] = 0.00000820413566106744f; //Vm millivolt
        sv[1] = 0.8789655121804799f;     //h dimensionless
    }
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) 
{

    // Get the mapping array
    uint32_t *mapping = NULL;
    if(extra_data) 
    {
        mapping = (uint32_t*)extra_data;
    }
    else 
    {
        print_to_stderr_and_file_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    uint32_t sv_id;

	int i;

    #pragma omp parallel for private(sv_id)
    for (i = 0; i < num_cells_to_solve; i++) 
    {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = (uint32_t )i;

        for (int j = 0; j < num_steps; ++j) 
        {
            if (mapping[i] == 0)
                solve_model_ode_cpu_fhn(dt, sv + (sv_id * NEQ_1), stim_currents[i]);
            else
                solve_model_ode_cpu_mitchell(dt, sv + (sv_id * NEQ_2), stim_currents[i]);
        }
    }
}

void solve_model_ode_cpu_fhn (real dt, real *sv, real stim_current)  
{

    real rY[NEQ_1], rDY[NEQ_1];

    for(int i = 0; i < NEQ_1; i++)
        rY[i] = sv[i];

    RHS_cpu_fhn(rY, rDY, stim_current);

    for(int i = 0; i < NEQ_1; i++)
        sv[i] = dt*rDY[i] + rY[i];
    
}

void RHS_cpu_fhn (const real *sv, real *rDY_, real stim_current) 
{

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

void solve_model_ode_cpu_mitchell (real dt, real *sv, real stim_current)
{
    real rY[NEQ_2], rDY[NEQ_2];

    for(int i = 0; i < NEQ_2; i++)
        rY[i] = sv[i];

    RHS_cpu_mitchell(rY, rDY, stim_current);

    for(int i = 0; i < NEQ_2; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu_mitchell(const real *sv, real *rDY_, real stim_current)
{
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