#include <assert.h>
#include <stdlib.h>
#include "ten_tusscher_tt3_mixed_endo_mid_epi.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;

}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using Mixed Ten & Tusscher 3 CPU model\n");

    uint32_t num_cells = solver->original_num_cells;

    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    OMP(parallel for)
        for(uint32_t i = 0; i < num_cells; i++) {

            real *sv = &solver->sv[i * NEQ];

            sv[0] = -86.2f;   // V;       millivolt
            sv[1] = 0.0f; //M
            sv[2] = 0.75; //H
            sv[3] = 0.75; //J
            sv[4] = 0.0f; //Xr1
            sv[5] = 0.0f; //Xs
            sv[6] = 1.0f; //S
            sv[7] = 1.0f; //F
            sv[8] = 1.0f; //F2
            sv[9] = 0.0; //D_INF
            sv[10] = 0.0; //R_INF
            sv[11] = 0.0; //Xr2_INF
        }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;
    real *transmurality = NULL;
    real *fibrosis = NULL;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    int num_extra_parameters = 8;
    real extra_par[num_extra_parameters];

    struct extra_data_for_tt3* extra_data_from_solver = (struct extra_data_for_tt3*)ode_solver->ode_extra_data;
    bool deallocate = false;

    if(ode_solver->ode_extra_data) {
        transmurality = extra_data_from_solver->transmurality;
        fibrosis = extra_data_from_solver->fibrosis;
        extra_par[0] = extra_data_from_solver->atpi;
        extra_par[1] = extra_data_from_solver->Ko;
        extra_par[2] = extra_data_from_solver->Ki;
        extra_par[3] = extra_data_from_solver->Vm_modifier;
        extra_par[4] = extra_data_from_solver->GNa_multiplicator;
        extra_par[5] = extra_data_from_solver->GCaL_multiplicator;
        extra_par[6] = extra_data_from_solver->INaCa_multiplicator;
        extra_par[7] = extra_data_from_solver->Ikatp_multiplicator;
    }
    else {
        // Default values for a healthy cell ///////////
        extra_par[0] = 6.8f;
        extra_par[1] = 5.4f;
        extra_par[2] = 138.3f;
        extra_par[3] = 0.0;
        extra_par[4] = 1.0f;
        extra_par[5] = 1.0f;
        extra_par[6] = 1.0f;
        extra_par[7] = 1.0f;

        fibrosis = (real*)malloc(num_cells_to_solve*sizeof(real));
        transmurality = (real*)malloc(num_cells_to_solve*sizeof(real));

        // Default case: all cells ENDO and Healthy
        for(uint64_t i = 0; i < num_cells_to_solve; i++) {
            fibrosis[i] = 1.0;          // Healthy
            transmurality[i] = 0.0;     // ENDO
        }

        deallocate = true;
    }

    int i;

    OMP(parallel for private(sv_id))
    for (i = 0; i < num_cells_to_solve; i++) {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], fibrosis[i], transmurality[i], extra_par);
        }
    }

    if(deallocate) {
        free(fibrosis);
        free(transmurality);
    }
}


void solve_model_ode_cpu(real dt, real *sv, real stim_current, real fibrosis, real transmurality, real *extra_parameters)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, fibrosis, transmurality, extra_parameters);

    //THIS MODEL USES THE Rush Larsen Method TO SOLVE THE EDOS
    sv[0] = dt*rDY[0] + rY[0];
    sv[1]  = rDY[1];
    sv[2]  = rDY[2];
    sv[3]  = rDY[3];
    sv[4]  = rDY[4];
    sv[5]  = rDY[5];
    sv[6]  = rDY[6];
    sv[7]  = rDY[7];
    sv[8]  = rDY[8];
    sv[9]  = rDY[9];
    sv[10]  = rDY[10];
    sv[11]  = rDY[11];
}


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real fibrosis, real transmurality, real const *extra_parameters) {

    //fibrosis = 0 means that the cell is fibrotic, 1 is not fibrotic. Anything between 0 and 1 means border zone
    //transmurality = 0 means cell is endocardium, 1 is mid-myocardium and 2 is epicardium

    //THIS IS THE STATE VECTOR THAT WE NEED TO SAVE IN THE STEADY STATE
    const real svolt    = sv[0];
    const real sm       = sv[1];
    const real sh       = sv[2];
    const real sj       = sv[3];
    const real sxr1     = sv[4];
    const real sxs      = sv[5];
    const real ss       = sv[6];
    const real sf       = sv[7];
    const real sf2      = sv[8];
    const real D_INF    = sv[9];
    const real R_INF    = sv[10];
    const real Xr2_INF  = sv[11];

    #include "ten_tusscher_tt3_mixed_endo_mid_epi.common.c"
}
