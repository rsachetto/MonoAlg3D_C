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

    log_info("Using Trovato_2019 CPU model\n");

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
    for(uint32_t i = 0; i < num_cells; i++) 
    {
        real *sv = &solver->sv[i * NEQ];
        
        // Steady-state 40 pulses (BCL=1000ms)
        sv[0] = -86.7099;
        sv[1] = 0.005431;
        sv[2] = 0.000104;
        sv[3] = 8.25533;
        sv[4] = 8.25502;
        sv[5] = 8.25503;
        sv[6] = 143.743;
        sv[7] = 143.744;
        sv[8] = 143.744;
        sv[9] = 4.4e-05;
        sv[10] = 0.000103;
        sv[11] = 1.26947;
        sv[12] = 1.25254;
        sv[13] = 1.27103;
        sv[14] = 1.1e-05;
        sv[15] = 0;
        sv[16] = 0.006303;
        sv[17] = 0.789469;
        sv[18] = 0.789392;
        sv[19] = 0.791301;
        sv[20] = 0.580955;
        sv[21] = 0.791719;
        sv[22] = 0.000241;
        sv[23] = 0.463851;
        sv[24] = 0.239936;
        sv[25] = 0.000272;
        sv[26] = 0.646362;
        sv[27] = 0.98999;
        sv[28] = 0;
        sv[29] = 1;
        sv[30] = 0.926919;
        sv[31] = 1;
        sv[32] = 1;
        sv[33] = 0.999976;
        sv[34] = 1;
        sv[35] = 1;
        sv[36] = 0.005885;
        sv[37] = 0.000303;
        sv[38] = 0.994251;
        sv[39] = 0.000367;
        sv[40] = 0.566131;
        sv[41] = 0.189842;
        sv[42] = 0.000222;
        sv[43] = 0.233515;
        sv[44] = 0.997077;
        sv[45] = 0.471259;

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

#include "../default_solvers.c"

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    const real V = sv[0];
    const real CaMKt = sv[1];
    const real cass = sv[2];
    const real nai = sv[3];
    const real nasl = sv[4];
    const real nass = sv[5];
    const real ki = sv[6];
    const real kss = sv[7];
    const real ksl = sv[8];
    const real cai = sv[9];
    const real casl = sv[10];
    const real cansr = sv[11];
    const real cajsr = sv[12];
    const real cacsr = sv[13];
    const real Jrel1 = sv[14];
    const real Jrel2 = sv[15];
    const real m = sv[16];
    const real hf = sv[17];
    const real hs = sv[18];
    const real j = sv[19];
    const real hsp = sv[20];
    const real jp = sv[21];
    const real mL = sv[22];
    const real hL = sv[23];
    const real hLp = sv[24];
    const real a = sv[25];
    const real i1 = sv[26];
    const real i2 = sv[27];
    const real d = sv[28];
    const real ff = sv[29];
    const real fs = sv[30];
    const real fcaf = sv[31];
    const real fcas = sv[32];
    const real jca = sv[33];
    const real ffp = sv[34];
    const real fcafp = sv[35];
    const real nca = sv[36];
    const real b = sv[37];
    const real g = sv[38];
    const real xrf = sv[39];
    const real xrs = sv[40];
    const real xs1 = sv[41];
    const real xs2 = sv[42];
    const real y = sv[43];
    const real xk1 = sv[44];
    const real u = sv[45]; 

    #include "trovato_2019_common.inc"
}
