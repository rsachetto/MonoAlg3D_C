#include "Maleckar2009.h"
#include <stdlib.h>

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using Maleckar2009 CPU model\n");

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

            sv[0] = -73.941851;
            sv[1] = 0.003325;
            sv[2] = 0.875262;
            sv[3] = 0.870692;
            sv[4] = 0.000014;
            sv[5] = 0.998578;
            sv[6] = 0.998561;
            sv[7] = 0.001098;
            sv[8] = 0.948202;
            sv[9] = 0.000371;
            sv[10] = 0.966869;
            sv[11] = 0.004661;
            sv[12] = 0.000054;
            sv[13] = 8.488527;
            sv[14] = 6.5e-5;
            sv[15] = 129.502075;
            sv[16] = 7.1e-5;
            sv[17] = 0.026604;
            sv[18] = 0.012843;
            sv[19] = 0.190077;
            sv[20] = 0.714719;
            sv[21] = 1.38222;
            sv[22] = 130.019282;
            sv[23] = 1.814418;
            sv[24] = 5.588239;
            sv[25] = 0.630471;
            sv[26] = 0.646226;
            sv[27] = 0.43071;
            sv[28] = 0.45453;
            sv[29] = 0.002665;

        }
}

void RHS_cpu(const real *sv, real *rDY, real stim_current, real dt) {

    real var_membrane__V = sv[0]; // Units: millivolt; Initial value: -73.941851
    real var_sodium_current_m_gate__m = sv[1]; // Units: dimensionless; Initial value: 0.003325
    real var_sodium_current_h1_gate__h1 = sv[2]; // Units: dimensionless; Initial value: 0.875262
    real var_sodium_current_h2_gate__h2 = sv[3]; // Units: dimensionless; Initial value: 0.870692
    real var_L_type_Ca_channel_d_L_gate__d_L = sv[4]; // Units: dimensionless; Initial value: 0.000014
    real var_L_type_Ca_channel_f_L1_gate__f_L1 = sv[5]; // Units: dimensionless; Initial value: 0.998578
    real var_L_type_Ca_channel_f_L2_gate__f_L2 = sv[6]; // Units: dimensionless; Initial value: 0.998561
    real var_Ca_independent_transient_outward_K_current_r_gate__r = sv[7]; // Units: dimensionless; Initial value: 0.001098
    real var_Ca_independent_transient_outward_K_current_s_gate__s = sv[8]; // Units: dimensionless; Initial value: 0.948202
    real var_ultra_rapid_K_current_aur_gate__a_ur = sv[9]; // Units: dimensionless; Initial value: 0.000371
    real var_ultra_rapid_K_current_iur_gate__i_ur = sv[10]; // Units: dimensionless; Initial value: 0.966869
    real var_delayed_rectifier_K_currents_n_gate__n = sv[11]; // Units: dimensionless; Initial value: 0.004661
    real var_delayed_rectifier_K_currents_pa_gate__pa = sv[12]; // Units: dimensionless; Initial value: 0.000054
    real var_intracellular_ion_concentrations__Na_i = sv[13]; // Units: millimolar; Initial value: 8.488527
    real var_intracellular_ion_concentrations__Ca_i = sv[14]; // Units: millimolar; Initial value: 6.5e-5
    real var_intracellular_ion_concentrations__K_i = sv[15]; // Units: millimolar; Initial value: 129.502075
    real var_intracellular_ion_concentrations__Ca_d = sv[16]; // Units: millimolar; Initial value: 7.1e-5
    real var_intracellular_Ca_buffering__O_C = sv[17]; // Units: dimensionless; Initial value: 0.026604
    real var_intracellular_Ca_buffering__O_TC = sv[18]; // Units: dimensionless; Initial value: 0.012843
    real var_intracellular_Ca_buffering__O_TMgC = sv[19]; // Units: dimensionless; Initial value: 0.190077
    real var_intracellular_Ca_buffering__O_TMgMg = sv[20]; // Units: dimensionless; Initial value: 0.714719
    real var_cleft_space_ion_concentrations__Na_c = sv[22]; // Units: millimolar; Initial value: 130.019282
    real var_cleft_space_ion_concentrations__Ca_c = sv[23]; // Units: millimolar; Initial value: 1.814418
    real var_cleft_space_ion_concentrations__K_c = sv[24]; // Units: millimolar; Initial value: 5.588239
    real var_Ca_handling_by_the_SR__Ca_rel = sv[25]; // Units: millimolar; Initial value: 0.630471
    real var_Ca_handling_by_the_SR__Ca_up = sv[26]; // Units: millimolar; Initial value: 0.646226
    real var_Ca_handling_by_the_SR__O_Calse = sv[27]; // Units: dimensionless; Initial value: 0.43071
    real var_Ca_handling_by_the_SR__F1 = sv[28]; // Units: dimensionless; Initial value: 0.45453
    real var_Ca_handling_by_the_SR__F2 = sv[29]; // Units: dimensionless; Initial value: 0.002665

    #include "Maleckar2009_common.inc.c"

}

#include "../default_solvers.c"
