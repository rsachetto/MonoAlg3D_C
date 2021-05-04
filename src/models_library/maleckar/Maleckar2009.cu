#include "Maleckar2009.h"
#include <stddef.h>
#include <stdint.h>

#define sv(i) *((real *)((char *)sv + pitch * (i)) + thread_id)

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {

	int thread_id = blockDim.x * blockIdx.x + threadIdx.x;

    if(thread_id < num_volumes) {

        sv(0) = -73.941851;
        sv(1) = 0.003325;
        sv(2) = 0.875262;
        sv(3) = 0.870692;
        sv(4) = 0.000014;
        sv(5) = 0.998578;
        sv(6) = 0.998561;
        sv(7) = 0.001098;
        sv(8) = 0.948202;
        sv(9) = 0.000371;
        sv(10) = 0.966869;
        sv(11) = 0.004661;
        sv(12) = 0.000054;
        sv(13) = 8.488527;
        sv(14) = 6.5e-5;
        sv(15) = 129.502075;
        sv(16) = 7.1e-5;
        sv(17) = 0.026604;
        sv(18) = 0.012843;
        sv(19) = 0.190077;
        sv(20) = 0.714719;
        sv(21) = 1.38222;
        sv(22) = 130.019282;
        sv(23) = 1.814418;
        sv(24) = 5.588239;
        sv(25) = 0.630471;
        sv(26) = 0.646226;
        sv(27) = 0.43071;
        sv(28) = 0.45453;
        sv(29) = 0.002665;

        if(use_adpt_dt) {
            sv(NEQ) = min_dt; // dt
            sv(NEQ+1) = 0.0;    // time_new
            sv(NEQ+2) = 0.0;    // previous dt
        }
    }
}

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int thread_id, real dt, size_t pitch, bool use_adpt_dt) {
    // State variables
    real var_membrane__V;                                          // Units: millivolt; Initial value: -73.941851
    real var_sodium_current_m_gate__m;                             // Units: dimensionless; Initial value: 0.003325
    real var_sodium_current_h1_gate__h1;                           // Units: dimensionless; Initial value: 0.875262
    real var_sodium_current_h2_gate__h2;                           // Units: dimensionless; Initial value: 0.870692
    real var_L_type_Ca_channel_d_L_gate__d_L;                      // Units: dimensionless; Initial value: 0.000014
    real var_L_type_Ca_channel_f_L1_gate__f_L1;                    // Units: dimensionless; Initial value: 0.998578
    real var_L_type_Ca_channel_f_L2_gate__f_L2;                    // Units: dimensionless; Initial value: 0.998561
    real var_Ca_independent_transient_outward_K_current_r_gate__r; // Units: dimensionless; Initial value: 0.001098
    real var_Ca_independent_transient_outward_K_current_s_gate__s; // Units: dimensionless; Initial value: 0.948202
    real var_ultra_rapid_K_current_aur_gate__a_ur;                 // Units: dimensionless; Initial value: 0.000371
    real var_ultra_rapid_K_current_iur_gate__i_ur;                 // Units: dimensionless; Initial value: 0.966869
    real var_delayed_rectifier_K_currents_n_gate__n;               // Units: dimensionless; Initial value: 0.004661
    real var_delayed_rectifier_K_currents_pa_gate__pa;             // Units: dimensionless; Initial value: 0.000054
    real var_intracellular_ion_concentrations__Na_i;               // Units: millimolar; Initial value: 8.488527
    real var_intracellular_ion_concentrations__Ca_i;               // Units: millimolar; Initial value: 6.5e-5
    real var_intracellular_ion_concentrations__K_i;                // Units: millimolar; Initial value: 129.502075
    real var_intracellular_ion_concentrations__Ca_d;               // Units: millimolar; Initial value: 7.1e-5
    real var_intracellular_Ca_buffering__O_C;                      // Units: dimensionless; Initial value: 0.026604
    real var_intracellular_Ca_buffering__O_TC;                     // Units: dimensionless; Initial value: 0.012843
    real var_intracellular_Ca_buffering__O_TMgC;                   // Units: dimensionless; Initial value: 0.190077
    real var_intracellular_Ca_buffering__O_TMgMg;                  // Units: dimensionless; Initial value: 0.714719
    real var_cleft_space_ion_concentrations__Na_c;                 // Units: millimolar; Initial value: 130.019282
    real var_cleft_space_ion_concentrations__Ca_c;                 // Units: millimolar; Initial value: 1.814418
    real var_cleft_space_ion_concentrations__K_c;                  // Units: millimolar; Initial value: 5.588239
    real var_Ca_handling_by_the_SR__Ca_rel;                        // Units: millimolar; Initial value: 0.630471
    real var_Ca_handling_by_the_SR__Ca_up;                         // Units: millimolar; Initial value: 0.646226
    real var_Ca_handling_by_the_SR__O_Calse;                       // Units: dimensionless; Initial value: 0.43071
    real var_Ca_handling_by_the_SR__F1;                            // Units: dimensionless; Initial value: 0.45453
    real var_Ca_handling_by_the_SR__F2;                            // Units: dimensionless; Initial value: 0.002665

    if(use_adpt_dt) {
        var_membrane__V = sv[0];                                          // Units: millivolt; Initial value: -73.941851
        var_sodium_current_m_gate__m = sv[1];                             // Units: dimensionless; Initial value: 0.003325
        var_sodium_current_h1_gate__h1 = sv[2];                           // Units: dimensionless; Initial value: 0.875262
        var_sodium_current_h2_gate__h2 = sv[3];                           // Units: dimensionless; Initial value: 0.870692
        var_L_type_Ca_channel_d_L_gate__d_L = sv[4];                      // Units: dimensionless; Initial value: 0.000014
        var_L_type_Ca_channel_f_L1_gate__f_L1 = sv[5];                    // Units: dimensionless; Initial value: 0.998578
        var_L_type_Ca_channel_f_L2_gate__f_L2 = sv[6];                    // Units: dimensionless; Initial value: 0.998561
        var_Ca_independent_transient_outward_K_current_r_gate__r = sv[7]; // Units: dimensionless; Initial value: 0.001098
        var_Ca_independent_transient_outward_K_current_s_gate__s = sv[8]; // Units: dimensionless; Initial value: 0.948202
        var_ultra_rapid_K_current_aur_gate__a_ur = sv[9];                 // Units: dimensionless; Initial value: 0.000371
        var_ultra_rapid_K_current_iur_gate__i_ur = sv[10];                // Units: dimensionless; Initial value: 0.966869
        var_delayed_rectifier_K_currents_n_gate__n = sv[11];              // Units: dimensionless; Initial value: 0.004661
        var_delayed_rectifier_K_currents_pa_gate__pa = sv[12];            // Units: dimensionless; Initial value: 0.000054
        var_intracellular_ion_concentrations__Na_i = sv[13];              // Units: millimolar; Initial value: 8.488527
        var_intracellular_ion_concentrations__Ca_i = sv[14];              // Units: millimolar; Initial value: 6.5e-5
        var_intracellular_ion_concentrations__K_i = sv[15];               // Units: millimolar; Initial value: 129.502075
        var_intracellular_ion_concentrations__Ca_d = sv[16];              // Units: millimolar; Initial value: 7.1e-5
        var_intracellular_Ca_buffering__O_C = sv[17];                     // Units: dimensionless; Initial value: 0.026604
        var_intracellular_Ca_buffering__O_TC = sv[18];                    // Units: dimensionless; Initial value: 0.012843
        var_intracellular_Ca_buffering__O_TMgC = sv[19];                  // Units: dimensionless; Initial value: 0.190077
        var_intracellular_Ca_buffering__O_TMgMg = sv[20];                 // Units: dimensionless; Initial value: 0.714719
        var_cleft_space_ion_concentrations__Na_c = sv[22];                // Units: millimolar; Initial value: 130.019282
        var_cleft_space_ion_concentrations__Ca_c = sv[23];                // Units: millimolar; Initial value: 1.814418
        var_cleft_space_ion_concentrations__K_c = sv[24];                 // Units: millimolar; Initial value: 5.588239
        var_Ca_handling_by_the_SR__Ca_rel = sv[25];                       // Units: millimolar; Initial value: 0.630471
        var_Ca_handling_by_the_SR__Ca_up = sv[26];                        // Units: millimolar; Initial value: 0.646226
        var_Ca_handling_by_the_SR__O_Calse = sv[27];                      // Units: dimensionless; Initial value: 0.43071
        var_Ca_handling_by_the_SR__F1 = sv[28];                           // Units: dimensionless; Initial value: 0.45453
        var_Ca_handling_by_the_SR__F2 = sv[29];                           // Units: dimensionless; Initial value: 0.002665

    } else {
        var_membrane__V = sv(0);                                          // Units: millivolt; Initial value: -73.941851
        var_sodium_current_m_gate__m = sv(1);                             // Units: dimensionless; Initial value: 0.003325
        var_sodium_current_h1_gate__h1 = sv(2);                           // Units: dimensionless; Initial value: 0.875262
        var_sodium_current_h2_gate__h2 = sv(3);                           // Units: dimensionless; Initial value: 0.870692
        var_L_type_Ca_channel_d_L_gate__d_L = sv(4);                      // Units: dimensionless; Initial value: 0.000014
        var_L_type_Ca_channel_f_L1_gate__f_L1 = sv(5);                    // Units: dimensionless; Initial value: 0.998578
        var_L_type_Ca_channel_f_L2_gate__f_L2 = sv(6);                    // Units: dimensionless; Initial value: 0.998561
        var_Ca_independent_transient_outward_K_current_r_gate__r = sv(7); // Units: dimensionless; Initial value: 0.001098
        var_Ca_independent_transient_outward_K_current_s_gate__s = sv(8); // Units: dimensionless; Initial value: 0.948202
        var_ultra_rapid_K_current_aur_gate__a_ur = sv(9);                 // Units: dimensionless; Initial value: 0.000371
        var_ultra_rapid_K_current_iur_gate__i_ur = sv(10);                // Units: dimensionless; Initial value: 0.966869
        var_delayed_rectifier_K_currents_n_gate__n = sv(11);              // Units: dimensionless; Initial value: 0.004661
        var_delayed_rectifier_K_currents_pa_gate__pa = sv(12);            // Units: dimensionless; Initial value: 0.000054
        var_intracellular_ion_concentrations__Na_i = sv(13);              // Units: millimolar; Initial value: 8.488527
        var_intracellular_ion_concentrations__Ca_i = sv(14);              // Units: millimolar; Initial value: 6.5e-5
        var_intracellular_ion_concentrations__K_i = sv(15);               // Units: millimolar; Initial value: 129.502075
        var_intracellular_ion_concentrations__Ca_d = sv(16);              // Units: millimolar; Initial value: 7.1e-5
        var_intracellular_Ca_buffering__O_C = sv(17);                     // Units: dimensionless; Initial value: 0.026604
        var_intracellular_Ca_buffering__O_TC = sv(18);                    // Units: dimensionless; Initial value: 0.012843
        var_intracellular_Ca_buffering__O_TMgC = sv(19);                  // Units: dimensionless; Initial value: 0.190077
        var_intracellular_Ca_buffering__O_TMgMg = sv(20);                 // Units: dimensionless; Initial value: 0.714719
        var_cleft_space_ion_concentrations__Na_c = sv(22);                // Units: millimolar; Initial value: 130.019282
        var_cleft_space_ion_concentrations__Ca_c = sv(23);                // Units: millimolar; Initial value: 1.814418
        var_cleft_space_ion_concentrations__K_c = sv(24);                 // Units: millimolar; Initial value: 5.588239
        var_Ca_handling_by_the_SR__Ca_rel = sv(25);                       // Units: millimolar; Initial value: 0.630471
        var_Ca_handling_by_the_SR__Ca_up = sv(26);                        // Units: millimolar; Initial value: 0.646226
        var_Ca_handling_by_the_SR__O_Calse = sv(27);                      // Units: dimensionless; Initial value: 0.43071
        var_Ca_handling_by_the_SR__F1 = sv(28);                           // Units: dimensionless; Initial value: 0.45453
        var_Ca_handling_by_the_SR__F2 = sv(29);                           // Units: dimensionless; Initial value: 0.002665
    }

	#include "Maleckar2009_common.inc.c"
}

//Include the default solver used by all models.
#include "../default_solvers.cu"
