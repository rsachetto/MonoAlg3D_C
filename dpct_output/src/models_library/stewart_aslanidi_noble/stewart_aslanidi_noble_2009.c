// TODO: Describe all the variable names

#include "stewart_aslanidi_noble_2009.h"
#include <stdio.h>

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu)
{

    log_info("Using Stewart-Aslanidi-Noble 2009 CPU model\n");

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

        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_info("Using Euler model to solve the ODEs\n");
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


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    real STATES[NEQ];
    for (uint32_t i = 0; i < NEQ; i++)
        STATES[i] = sv[i];

    #include "stewart_aslanidi_noble_2009_common.inc"
}

#include "../default_solvers.c"
