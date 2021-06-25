#include <assert.h>
#include <stdlib.h>
#include "ten_tusscher_2004_mixed_endo_mid_epi.h"


GET_CELL_MODEL_DATA(init_cell_model_data) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;

}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {

        real *sv = &solver->sv[i * NEQ];

        sv[0] = INITIAL_V; // V;       millivolt
        sv[1] = 0.f;       // M
        sv[2] = 0.75;      // H
        sv[3] = 0.75f;     // J
        sv[4] = 0.f;       // Xr1
        sv[5] = 1.f;       // Xr2
        sv[6] = 0.f;       // Xs
        sv[7] = 1.f;       // S
        sv[8] = 0.f;       // R
        sv[9] = 0.f;       // D
        sv[10] = 1.f;      // F
        sv[11] = 1.f;      // FCa
        sv[12] = 1.f;      // G
        sv[13] = 0.0002;   // Cai
        sv[14] = 0.2f;     // CaSR
        sv[15] = 11.6f;    // Nai
        sv[16] = 138.3f;   // Ki
    }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;
    int i;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;


    // Get the mapping array
    uint32_t *mapping = NULL;
    if(ode_solver->ode_extra_data) {
        mapping = (uint32_t*)ode_solver->ode_extra_data;
    }
    else {
        log_error_and_exit("You need to specify a mask function when using this mixed model!\n");
    }


    OMP(parallel for private(sv_id))
    for (i = 0; i < num_cells_to_solve; i++) {

        if(ode_solver->cells_to_solve)
            sv_id = ode_solver->cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < ode_solver->num_steps; ++j) {
            solve_model_ode_cpu(ode_solver->min_dt, ode_solver->sv + (sv_id * NEQ), stim_currents[i], mapping[i]);
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current, int mapping)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, mapping);

    for(int i = 0; i < NEQ; i++)
        sv[i] = rDY[i];
}


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, int mapping) {

    // State variables
    real svolt = sv[0];
    real sm    = sv[1];
    real sh    = sv[2];
    real sj    = sv[3];
    real sxr1  = sv[4];
    real sxr2  = sv[5];
    real sxs   = sv[6];
    real ss    = sv[7];
    real sr    = sv[8];
    real sd    = sv[9];
    real sf    = sv[10];
    real sfca  = sv[11];
    real sg    = sv[12];
    real Cai   = sv[13];
    real CaSR  = sv[14];
    real Nai   = sv[15];
    real Ki    = sv[16];

    #include "ten_tusscher_2004_mixed_endo_mid_epi.common.c"

}
