#include "model_common.h"
#include <assert.h>
#include <stdlib.h>
#include "ten_tusscher_3_RS.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;

}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

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
    sv[11] = 0.0; //Xr2_INF}
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    uint32_t sv_id;
    real *fibrosis;

    // Default values for a healthy cell ///////////
    real atpi = 6.8f;
    real Ko = 5.4f;
    real Ki_mult = 1.0f;
    real acidosis = 0.0;
    real K1_mult = 1.0f;
    ////////////////////////////////////

    if(extra_data) {
        atpi = ((real*)extra_data)[0]; //value
        Ko = ((real*)extra_data)[1]; //value
        Ki_mult = ((real*)extra_data)[2]; //value
        K1_mult = ((real*)extra_data)[3]; //value
        acidosis = ((real*)extra_data)[4]; //value
        fibrosis = ((real*)extra_data) + 5; //pointer
    }
    else {
        fibrosis = calloc(num_cells_to_solve, sizeof(real));
    }

	int i;

    #pragma omp parallel for private(sv_id)
    for (i = 0; i < num_cells_to_solve; i++) {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], fibrosis[i], atpi, Ko, Ki_mult, K1_mult, acidosis);
        }
    }

    if(extra_data == NULL) free(fibrosis);
}


void solve_model_ode_cpu(real dt, real *sv, real stim_current, real fibrosis, real atpi, real Ko, real Ki_mult,
                         real K1_mult, real acidosis)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, fibrosis, atpi, Ko, Ki_mult, K1_mult, acidosis);

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


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real fibrosis, real atpi, real Ko,
             real Ki_multiplicator, real K1_multiplicator, real acidosis) {

    // State variables
    const real svolt = sv[0];      // Membrane variable

    real svolt_acid = svolt;

    if( (fibrosis == 0.0f) && (acidosis == 1.0f) ) {
        //These values are from In Electrophysiologic effects of acute myocardial ischemia: a theoretical
        //study of altered cell excitability and action potential duration
        svolt_acid = svolt - 3.4f;
    }

    const real sm   = sv[1];
    const real sh   = sv[2];
    const real sj   = sv[3];
    const real sxr1 = sv[4];
    const real sxs  = sv[5];
    const real ss   = sv[6];
    const real  sf   = sv[7];
    const real sf2  = sv[8];

    const real D_INF  = sv[9];
    const real Xr2_INF  = sv[10];
    const real R_INF  = sv[11];

    #include "ten_tusscher_3_common.inc"

}
