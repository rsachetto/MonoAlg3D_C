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


	char *cell_type;
#ifdef ENDO
	cell_type = strdup("ENDO");
#endif

#ifdef EPI
	cell_type = strdup("EPI");
#endif

#ifdef MCELL
	cell_type = strdup("MCELL");
#endif

	log_to_stdout_and_file("Using ten Tusscher 3 %s CPU model\n", cell_type);
	free(cell_type);

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
    real *fibrosis;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    // Default values for a healthy cell ///////////
    real atpi = 6.8f;
    real Ko = 5.4f;
    real Ki = 138.3f;
    real Vm_change = 0.0;
    real GNa_multiplicator = 1.0f;
    real GCa_multiplicator = 1.0f;
    ////////////////////////////////////

    int num_extra_parameters = 6;
    size_t extra_parameters_size = num_extra_parameters*sizeof(real);

    if(ode_solver->ode_extra_data) {
        fibrosis = ((real*)ode_solver->ode_extra_data) + num_extra_parameters; //pointer
    }
    else {
        ode_solver->ode_extra_data = malloc(extra_parameters_size);
        ((real*)ode_solver->ode_extra_data)[0] = atpi;
        ((real*)ode_solver->ode_extra_data)[1] = Ko;
        ((real*)ode_solver->ode_extra_data)[2] = Ki;
        ((real*)ode_solver->ode_extra_data)[3] = Vm_change;
        ((real*)ode_solver->ode_extra_data)[4] = GNa_multiplicator;
        ((real*)ode_solver->ode_extra_data)[5] = GCa_multiplicator;

        fibrosis = calloc(num_cells_to_solve, sizeof(real));
    }

    int i;

    OMP(parallel for private(sv_id))
    for (i = 0; i < num_cells_to_solve; i++) {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], fibrosis[i], ode_solver->ode_extra_data);

        }
    }

    if(ode_solver->ode_extra_data == NULL) free(fibrosis);
}


void solve_model_ode_cpu(real dt, real *sv, real stim_current, real fibrosis, real *extra_parameters)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, fibrosis, extra_parameters);

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


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real fibrosis, real *extra_parameters) {

    //fibrosis = 0 means that the cell is fibrotic, 1 is not fibrotic. Anything between 0 and 1 means border zone
        
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

    #include "ten_tusscher_3_RS_common.inc"
}
