#include <assert.h>
#include <stdlib.h>
#include "tt2004_myo_epi.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    static bool first_call = true;

    if(first_call) {
        log_to_stdout_and_file("Using mixed version of TT3 (MCELL + EPI) CPU model\n");
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
    solver->sv = (real*)malloc(NEQ*num_volumes*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_volumes; i++) {

        // MCELL
        if (mapping[i] == 0) {
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
        // EPI
        else if (mapping[i] == 1) {
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

            // MCELL
            if (mapping[i] == 0) {
                solve_model_ode_cpu_myo(dt, sv + (sv_id * NEQ), stim_currents[i]);
            }
            // EPI
            else if (mapping[i] == 1) {
                solve_model_ode_cpu_epi(dt, sv + (sv_id * NEQ), stim_currents[i]);
            }
            else {
                log_to_stderr_and_file_and_exit("[ERR] Invalid mapping value '%u'! It should be 0(MCELL) or 1(EPI)!\n",mapping[i]);
            }
            
        }
    }

}

void solve_model_ode_cpu_myo (real dt, real *sv, real stim_current)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu_myo(rY, rDY, stim_current, dt);

    for(int i = 0; i < NEQ; i++)
        sv[i] = rDY[i];
    
}

void RHS_cpu_myo (const real *sv, real *rDY_, real stim_current, real dt) {

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

    // Specific MCELL parameters
    real Gks = 0.062;
    real Gto = 0.294;

    real R_INF=1./(1.+exp((20-svolt)/6.));
    real S_INF=1./(1.+exp((svolt+20)/5.));
    real TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    real TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;

    #include "tt2004_common.inc"

}

void solve_model_ode_cpu_epi (real dt, real *sv, real stim_current) {
    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu_epi(rY, rDY, stim_current, dt);

    for(int i = 0; i < NEQ; i++)
        sv[i] = rDY[i];
}

void RHS_cpu_epi(const real *sv, real *rDY_, real stim_current, real dt) {

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

    // Specific EPI parameters
    real Gks = 0.245;
    real Gto = 0.294;

    real R_INF=1./(1.+exp((20-svolt)/6.));
    real S_INF=1./(1.+exp((svolt+20)/5.));
    real TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    real TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;

    #include "tt2004_common.inc"

}