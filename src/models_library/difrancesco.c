// A Model of Cardiac Electrical Activity Incorporating Ionic Pumps and Concentration Changes 
// Reference: https://models.physiomeproject.org/exposure/91d93b61d7da56b6baf1f0c4d88ecd77/difrancesco_noble_1985.cellml

#include "difrancesco.h"

real CONSTANTS[50];

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_to_stdout_and_file("Using DiFrancesco & Noble 1985 CPU model\n");

    uint32_t num_volumes = solver->original_num_cells;
	
	solver->sv = (real*)malloc(NEQ*num_volumes*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_volumes; i++) {
        real *sv = &solver->sv[i * NEQ];
        // Normal
        sv[0] = -87;    // V millivolt
        sv[1] = 4;      // Kc millimolar
        sv[2] = 140;    // Ki millimolar
        sv[3] = 8;      // Nai millimolar
        sv[4] = 0.2;    // y dimensionless
        sv[5] = 0.01;   // x dimensionless
        sv[6] = 5e-5;   // Cai millimolar
        sv[7] = 1;      // s dimensionless
        sv[8] = 0.01;   // m dimensionless
        sv[9] = 0.8;    // h dimensionless
        sv[10] = 0.005; // d dimensionless
        sv[11] = 1;     // f dimensionless
        sv[12] = 1;     // f2 dimensionless
        sv[13] = 2;     // Ca_up millimolar
        sv[14] = 1;     // Ca_rel millimolar
        sv[15] = 1;     // p dimensionless
    }

    CONSTANTS[0] = 8314.472;
    CONSTANTS[1] = 310;
    CONSTANTS[2] = 96485.3415;
    CONSTANTS[3] = 0.075;
    CONSTANTS[4] = 0;
    CONSTANTS[5] = 3;
    CONSTANTS[6] = 3;
    CONSTANTS[7] = 45;
    CONSTANTS[8] = 140;
    CONSTANTS[9] = 1e-5;
    CONSTANTS[10] = 180;
    CONSTANTS[11] = 920;
    CONSTANTS[12] = 210;
    CONSTANTS[13] = 10;
    CONSTANTS[14] = 0.0005;
    CONSTANTS[15] = 0.28;
    CONSTANTS[16] = 0.18;
    CONSTANTS[17] = 0.02;
    CONSTANTS[18] = 2;
    CONSTANTS[19] = 125;
    CONSTANTS[20] = 1;
    CONSTANTS[21] = 40;
    CONSTANTS[22] = 3;
    CONSTANTS[23] = 0.02;
    CONSTANTS[24] = 0.001;
    CONSTANTS[25] = 0.5;
    CONSTANTS[26] = 750;
    CONSTANTS[27] = 1e-5;
    CONSTANTS[28] = 15;
    CONSTANTS[29] = 0.0001;
    CONSTANTS[30] = 0.0001;
    CONSTANTS[31] = 5;
    CONSTANTS[32] = 0.001;
    CONSTANTS[33] = 0.05;
    CONSTANTS[34] = 2;
    CONSTANTS[35] = 0.1;
    CONSTANTS[36] = 5;
    CONSTANTS[37] = 0.001;
    CONSTANTS[38] = 0.025;
    CONSTANTS[39] = 2;
    CONSTANTS[40] = 0.05;
    CONSTANTS[41] = 2;
    CONSTANTS[42] = 0.00157;
    CONSTANTS[43] = 4;
    CONSTANTS[44] = 0.7;
    CONSTANTS[45] = ( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2];
    CONSTANTS[46] =  3.14159*pow(CONSTANTS[33], 2.00000)*CONSTANTS[34];
    CONSTANTS[47] =  CONSTANTS[46]*(1.00000 - CONSTANTS[35]);
    CONSTANTS[48] =  CONSTANTS[47]*0.0500000;
    CONSTANTS[49] = CONSTANTS[47]*0.0200000;
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    OMP(parallel for private(sv_id))
    for (uint32_t i = 0; i < num_cells_to_solve; i++)
    {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) 
        {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  
{

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current);

    // Solve the model using Forward Euler
    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current) {

    CONSTANTS[4] = stim_current;

    real STATES[16]; 
    STATES[0] = sv[0];
    STATES[1] = sv[1];
    STATES[2] = sv[2];
    STATES[3] = sv[3];
    STATES[4] = sv[4];
    STATES[5] = sv[5];
    STATES[6] = sv[6];
    STATES[7] = sv[7];
    STATES[8] = sv[8];
    STATES[9] = sv[9];
    STATES[10] = sv[10];
    STATES[11] = sv[11];
    STATES[12] = sv[12];
    STATES[13] = sv[13];
    STATES[14] = sv[14];
    STATES[15] = sv[15];


    double ALGEBRAIC[46];
    ALGEBRAIC[8] = ( STATES[6]*CONSTANTS[31])/CONSTANTS[32];
    ALGEBRAIC[2] = ( 0.500000*exp( 0.0826000*(STATES[0]+50.0000)))/(1.00000+exp( 0.0570000*(STATES[0]+50.0000)));
    ALGEBRAIC[12] = ( 1.30000*exp( - 0.0600000*(STATES[0]+20.0000)))/(1.00000+exp( - 0.0400000*(STATES[0]+20.0000)));
    ALGEBRAIC[3] =  0.0330000*exp(- STATES[0]/17.0000);
    ALGEBRAIC[13] = 33.0000/(1.00000+exp(- (STATES[0]+10.0000)/8.00000));
    ALGEBRAIC[5] =  20.0000*exp( - 0.125000*(STATES[0]+75.0000));
    ALGEBRAIC[15] = 2000.00/( 320.000*exp( - 0.100000*(STATES[0]+75.0000))+1.00000);
    ALGEBRAIC[9] = ( 0.625000*(STATES[0]+34.0000))/(exp((STATES[0]+34.0000)/4.00000) - 1.00000);
    ALGEBRAIC[18] = 5.00000/(1.00000+exp(( - 1.00000*(STATES[0]+34.0000))/4.00000));
    ALGEBRAIC[1] =  0.0500000*exp( - 0.0670000*((STATES[0]+52.0000) - 10.0000));
    ALGEBRAIC[11] = (STATES[0]+52.0000) - 10.0000;
    ALGEBRAIC[19] = (fabs(ALGEBRAIC[11])<CONSTANTS[9] ? 2.50000 : ( 1.00000*ALGEBRAIC[11])/(1.00000 - exp( - 0.200000*ALGEBRAIC[11])));
    ALGEBRAIC[4] = STATES[0]+41.0000;
    ALGEBRAIC[14] = (fabs(ALGEBRAIC[4])<CONSTANTS[27] ? 2000.00 : ( 200.000*ALGEBRAIC[4])/(1.00000 - exp( - 0.100000*ALGEBRAIC[4])));
    ALGEBRAIC[20] =  8000.00*exp( - 0.0560000*(STATES[0]+66.0000));
    ALGEBRAIC[6] = (STATES[0]+24.0000) - 5.00000;
    ALGEBRAIC[16] = (fabs(ALGEBRAIC[6])<CONSTANTS[29] ? 120.000 : ( 30.0000*ALGEBRAIC[6])/(1.00000 - exp(( - 1.00000*ALGEBRAIC[6])/4.00000)));
    ALGEBRAIC[21] = (fabs(ALGEBRAIC[6])<CONSTANTS[29] ? 120.000 : ( 12.0000*ALGEBRAIC[6])/(exp(ALGEBRAIC[6]/10.0000) - 1.00000));
    ALGEBRAIC[7] = STATES[0]+34.0000;
    ALGEBRAIC[17] = (fabs(ALGEBRAIC[7])<CONSTANTS[30] ? 25.0000 : ( 6.25000*ALGEBRAIC[7])/(exp(ALGEBRAIC[7]/4.00000) - 1.00000));
    ALGEBRAIC[22] = 50.0000/(1.00000+exp(( - 1.00000*(STATES[0]+34.0000))/4.00000));
    ALGEBRAIC[0] =  CONSTANTS[45]*log(CONSTANTS[8]/STATES[3]);
    ALGEBRAIC[30] =  CONSTANTS[16]*(STATES[0] - ALGEBRAIC[0]);
    ALGEBRAIC[33] = ( (( CONSTANTS[19]*STATES[1])/(CONSTANTS[20]+STATES[1]))*STATES[3])/(CONSTANTS[21]+STATES[3]);
    ALGEBRAIC[34] = ( CONSTANTS[23]*( exp(( CONSTANTS[25]*(CONSTANTS[22] - 2.00000)*STATES[0])/CONSTANTS[45])*pow(STATES[3], CONSTANTS[22])*CONSTANTS[18] -  exp(( (CONSTANTS[25] - 1.00000)*(CONSTANTS[22] - 2.00000)*STATES[0])/CONSTANTS[45])*pow(CONSTANTS[8], CONSTANTS[22])*STATES[6]))/( (1.00000+ CONSTANTS[24]*( STATES[6]*pow(CONSTANTS[8], CONSTANTS[22])+ CONSTANTS[18]*pow(STATES[3], CONSTANTS[22])))*(1.00000+STATES[6]/0.00690000));
    ALGEBRAIC[35] =  CONSTANTS[45]*log((CONSTANTS[8]+ 0.120000*STATES[1])/(STATES[3]+ 0.120000*STATES[2]));
    ALGEBRAIC[36] =  CONSTANTS[26]*pow(STATES[8], 3.00000)*STATES[9]*(STATES[0] - ALGEBRAIC[35]);
    ALGEBRAIC[23] =  (( STATES[4]*STATES[1])/(STATES[1]+CONSTANTS[7]))*CONSTANTS[5]*(STATES[0] - ALGEBRAIC[0]);
    ALGEBRAIC[40] =  (( 0.0100000*CONSTANTS[28]*(STATES[0] - 50.0000))/( CONSTANTS[45]*(1.00000 - exp(( - 1.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))))*( STATES[3]*exp(50.0000/CONSTANTS[45]) -  CONSTANTS[8]*exp(( - 1.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))*STATES[10]*STATES[11]*STATES[12];
    ALGEBRAIC[39] =  (( 2.00000*1.00000*CONSTANTS[47]*CONSTANTS[2])/( 1.00000*CONSTANTS[38]*CONSTANTS[36]))*STATES[6]*(CONSTANTS[36] - STATES[13]);
    ALGEBRAIC[41] =  (( 2.00000*1.00000*CONSTANTS[49]*CONSTANTS[2])/( 1.00000*CONSTANTS[39]))*STATES[15]*(STATES[13] - STATES[14]);
    ALGEBRAIC[26] = ( CONSTANTS[10]*(STATES[2] -  STATES[1]*exp(- STATES[0]/CONSTANTS[45])))/140.000;
    ALGEBRAIC[27] =  STATES[5]*ALGEBRAIC[26];
    ALGEBRAIC[10] =  CONSTANTS[45]*log(STATES[1]/STATES[2]);
    ALGEBRAIC[28] = ( (( CONSTANTS[11]*STATES[1])/(STATES[1]+CONSTANTS[12]))*(STATES[0] - ALGEBRAIC[10]))/(1.00000+exp(( ((STATES[0]+10.0000) - ALGEBRAIC[10])*2.00000)/CONSTANTS[45]));
    ALGEBRAIC[29] =  (( (( STATES[7]*CONSTANTS[15]*(0.200000+STATES[1]/(CONSTANTS[13]+STATES[1]))*STATES[6])/(CONSTANTS[14]+STATES[6]))*(STATES[0]+10.0000))/(1.00000 - exp( - 0.200000*(STATES[0]+10.0000))))*( STATES[2]*exp(( 0.500000*STATES[0])/CONSTANTS[45]) -  STATES[1]*exp(( - 0.500000*STATES[0])/CONSTANTS[45]));
    ALGEBRAIC[24] =  (( STATES[4]*STATES[1])/(STATES[1]+CONSTANTS[7]))*CONSTANTS[6]*(STATES[0] - ALGEBRAIC[10]);
    ALGEBRAIC[38] =  (( 0.0100000*CONSTANTS[28]*(STATES[0] - 50.0000))/( CONSTANTS[45]*(1.00000 - exp(( - 1.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))))*( STATES[2]*exp(50.0000/CONSTANTS[45]) -  STATES[1]*exp(( - 1.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))*STATES[10]*STATES[11]*STATES[12];
    ALGEBRAIC[42] = (ALGEBRAIC[28]+ALGEBRAIC[27]+ALGEBRAIC[24]+ALGEBRAIC[38]+ALGEBRAIC[29]) -  2.00000*ALGEBRAIC[33];
    ALGEBRAIC[25] = ALGEBRAIC[23]+ALGEBRAIC[24];
    ALGEBRAIC[31] =  0.500000*CONSTANTS[45]*log(CONSTANTS[18]/STATES[6]);
    ALGEBRAIC[32] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[31]);
    ALGEBRAIC[37] =  (( 4.00000*CONSTANTS[28]*(STATES[0] - 50.0000))/( CONSTANTS[45]*(1.00000 - exp(( - 1.00000*(STATES[0] - 50.0000)*2.00000)/CONSTANTS[45]))))*( STATES[6]*exp(100.000/CONSTANTS[45]) -  CONSTANTS[18]*exp(( - 2.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))*STATES[10]*STATES[11]*STATES[12];
    ALGEBRAIC[43] = ALGEBRAIC[37]+ALGEBRAIC[38]+ALGEBRAIC[40];
    ALGEBRAIC[44] = ( (( 2.00000*1.00000*CONSTANTS[49]*CONSTANTS[2])/( 1.00000*CONSTANTS[40]))*STATES[14]*pow(STATES[6], CONSTANTS[41]))/(pow(STATES[6], CONSTANTS[41])+pow(CONSTANTS[37], CONSTANTS[41]));

    double RATES[16];
    RATES[0] = (- (ALGEBRAIC[25]+ALGEBRAIC[27]+ALGEBRAIC[28]+ALGEBRAIC[29]+ALGEBRAIC[30]+ALGEBRAIC[32]+ALGEBRAIC[33]+ALGEBRAIC[34]+ALGEBRAIC[36]+ALGEBRAIC[43]+CONSTANTS[4])/CONSTANTS[3]) * 1.0E-03;
    RATES[1] =  (- CONSTANTS[44]*(STATES[1] - CONSTANTS[43])+( 1.00000*ALGEBRAIC[42])/( 1.00000*CONSTANTS[42]*CONSTANTS[2])) * 1.0E-03;
    RATES[2] = ( ( - 1.00000*ALGEBRAIC[42])/( 1.00000*CONSTANTS[47]*CONSTANTS[2]) ) * 1.0E-03;
    RATES[3] = ( ( - 1.00000*(ALGEBRAIC[36]+ALGEBRAIC[30]+ALGEBRAIC[23]+ALGEBRAIC[40]+ ALGEBRAIC[33]*3.00000+( ALGEBRAIC[34]*CONSTANTS[22])/(CONSTANTS[22] - 2.00000)))/( 1.00000*CONSTANTS[47]*CONSTANTS[2])) * 1.0E-03;
    RATES[4] =  (ALGEBRAIC[1]*(1.00000 - STATES[4]) -  ALGEBRAIC[19]*STATES[4]) * 1.0E-03;
    RATES[5] =  (ALGEBRAIC[2]*(1.00000 - STATES[5]) -  ALGEBRAIC[12]*STATES[5]) * 1.0E-03;
    RATES[6] = (( - 1.00000*((((ALGEBRAIC[37]+ALGEBRAIC[32]) - ( 2.00000*ALGEBRAIC[34])/(CONSTANTS[22] - 2.00000)) - ALGEBRAIC[44])+ALGEBRAIC[39]))/( 2.00000*1.00000*CONSTANTS[47]*CONSTANTS[2])) * 1.0E-03;
    RATES[7] =  (ALGEBRAIC[3]*(1.00000 - STATES[7]) -  ALGEBRAIC[13]*STATES[7]) * 1.0E-03;
    RATES[8] =  (ALGEBRAIC[14]*(1.00000 - STATES[8]) -  ALGEBRAIC[20]*STATES[8]) * 1.0E-03;
    RATES[9] =  (ALGEBRAIC[5]*(1.00000 - STATES[9]) -  ALGEBRAIC[15]*STATES[9]) * 1.0E-03;
    RATES[10] =  (ALGEBRAIC[16]*(1.00000 - STATES[10]) -  ALGEBRAIC[21]*STATES[10]) * 1.0E-03;
    RATES[11] =  (ALGEBRAIC[17]*(1.00000 - STATES[11]) -  ALGEBRAIC[22]*STATES[11]) * 1.0E-03;
    RATES[12] = (CONSTANTS[31] -  STATES[12]*(CONSTANTS[31]+ALGEBRAIC[8])) * 1.0E-03;
    RATES[13] = (( 1.00000*(ALGEBRAIC[39] - ALGEBRAIC[41]))/( 2.00000*1.00000*CONSTANTS[48]*CONSTANTS[2])) * 1.0E-03;
    RATES[14] = (( 1.00000*(ALGEBRAIC[41] - ALGEBRAIC[44]))/( 2.00000*1.00000*CONSTANTS[49]*CONSTANTS[2])) * 1.0E-03;
    RATES[15] = (ALGEBRAIC[9]*(1.00000 - STATES[15]) - ALGEBRAIC[18]*STATES[15]) * 1.0E-03;


    rDY_[0] = RATES[0];
    rDY_[1] = RATES[1];
    rDY_[2] = RATES[2];
    rDY_[3] = RATES[3];
    rDY_[4] = RATES[4];
    rDY_[5] = RATES[5];
    rDY_[6] = RATES[6];
    rDY_[7] = RATES[7];
    rDY_[8] = RATES[8];
    rDY_[9] = RATES[9];
    rDY_[10] = RATES[10];
    rDY_[11] = RATES[11];
    rDY_[12] = RATES[12];
    rDY_[13] = RATES[13];
    rDY_[14] = RATES[14];
    rDY_[15] = RATES[15];

}

/*
 * VOI is time in component environment (second).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (joule_per_kilomole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_mole).
 * CONSTANTS[45] is RTONF in component membrane (millivolt).
 * CONSTANTS[3] is C in component membrane (microF).
 * CONSTANTS[4] is i_pulse in component membrane (nanoA).
 * ALGEBRAIC[25] is i_f in component hyperpolarising_activated_current (nanoA).
 * ALGEBRAIC[27] is i_K in component time_dependent_potassium_current (nanoA).
 * ALGEBRAIC[28] is i_K1 in component time_independent_potassium_current (nanoA).
 * ALGEBRAIC[29] is i_to in component transient_outward_current (nanoA).
 * ALGEBRAIC[30] is i_Na_b in component sodium_background_current (nanoA).
 * ALGEBRAIC[32] is i_Ca_b in component calcium_background_current (nanoA).
 * ALGEBRAIC[33] is i_p in component sodium_potassium_pump (nanoA).
 * ALGEBRAIC[34] is i_NaCa in component Na_Ca_exchanger (nanoA).
 * ALGEBRAIC[36] is i_Na in component fast_sodium_current (nanoA).
 * ALGEBRAIC[43] is i_si in component second_inward_current (nanoA).
 * ALGEBRAIC[23] is i_fNa in component hyperpolarising_activated_current (nanoA).
 * ALGEBRAIC[0] is E_Na in component hyperpolarising_activated_current (millivolt).
 * ALGEBRAIC[10] is E_K in component hyperpolarising_activated_current (millivolt).
 * ALGEBRAIC[24] is i_fK in component hyperpolarising_activated_current (nanoA).
 * CONSTANTS[5] is g_f_Na in component hyperpolarising_activated_current (microS).
 * CONSTANTS[6] is g_f_K in component hyperpolarising_activated_current (microS).
 * CONSTANTS[7] is Km_f in component hyperpolarising_activated_current (millimolar).
 * STATES[1] is Kc in component extracellular_potassium_concentration (millimolar).
 * STATES[2] is Ki in component intracellular_potassium_concentration (millimolar).
 * STATES[3] is Nai in component intracellular_sodium_concentration (millimolar).
 * CONSTANTS[8] is Nao in component extracellular_sodium_concentration (millimolar).
 * STATES[4] is y in component hyperpolarising_activated_current_y_gate (dimensionless).
 * ALGEBRAIC[1] is alpha_y in component hyperpolarising_activated_current_y_gate (per_second).
 * ALGEBRAIC[19] is beta_y in component hyperpolarising_activated_current_y_gate (per_second).
 * CONSTANTS[9] is delta_y in component hyperpolarising_activated_current_y_gate (millivolt).
 * ALGEBRAIC[11] is E0_y in component hyperpolarising_activated_current_y_gate (millivolt).
 * ALGEBRAIC[26] is I_K in component time_dependent_potassium_current (nanoA).
 * CONSTANTS[10] is i_K_max in component time_dependent_potassium_current (nanoA).
 * STATES[5] is x in component time_dependent_potassium_current_x_gate (dimensionless).
 * ALGEBRAIC[2] is alpha_x in component time_dependent_potassium_current_x_gate (per_second).
 * ALGEBRAIC[12] is beta_x in component time_dependent_potassium_current_x_gate (per_second).
 * CONSTANTS[11] is g_K1 in component time_independent_potassium_current (microS).
 * CONSTANTS[12] is Km_K1 in component time_independent_potassium_current (millimolar).
 * CONSTANTS[13] is Km_to in component transient_outward_current (millimolar).
 * CONSTANTS[14] is Km_Ca in component transient_outward_current (millimolar).
 * CONSTANTS[15] is g_to in component transient_outward_current (microS_per_millimolar).
 * STATES[6] is Cai in component intracellular_calcium_concentration (millimolar).
 * STATES[7] is s in component transient_outward_current_s_gate (dimensionless).
 * ALGEBRAIC[3] is alpha_s in component transient_outward_current_s_gate (per_second).
 * ALGEBRAIC[13] is beta_s in component transient_outward_current_s_gate (per_second).
 * CONSTANTS[16] is g_Nab in component sodium_background_current (microS).
 * ALGEBRAIC[31] is E_Ca in component calcium_background_current (millivolt).
 * CONSTANTS[17] is g_Cab in component calcium_background_current (microS).
 * CONSTANTS[18] is Cao in component extracellular_calcium_concentration (millimolar).
 * CONSTANTS[19] is I_p in component sodium_potassium_pump (nanoA).
 * CONSTANTS[20] is K_mK in component sodium_potassium_pump (millimolar).
 * CONSTANTS[21] is K_mNa in component sodium_potassium_pump (millimolar).
 * CONSTANTS[22] is n_NaCa in component Na_Ca_exchanger (dimensionless).
 * CONSTANTS[23] is K_NaCa in component Na_Ca_exchanger (nanoA).
 * CONSTANTS[24] is d_NaCa in component Na_Ca_exchanger (dimensionless).
 * CONSTANTS[25] is gamma in component Na_Ca_exchanger (dimensionless).
 * CONSTANTS[26] is g_Na in component fast_sodium_current (microS).
 * ALGEBRAIC[35] is E_mh in component fast_sodium_current (millivolt).
 * STATES[8] is m in component fast_sodium_current_m_gate (dimensionless).
 * STATES[9] is h in component fast_sodium_current_h_gate (dimensionless).
 * ALGEBRAIC[14] is alpha_m in component fast_sodium_current_m_gate (per_second).
 * ALGEBRAIC[20] is beta_m in component fast_sodium_current_m_gate (per_second).
 * CONSTANTS[27] is delta_m in component fast_sodium_current_m_gate (millivolt).
 * ALGEBRAIC[4] is E0_m in component fast_sodium_current_m_gate (millivolt).
 * ALGEBRAIC[5] is alpha_h in component fast_sodium_current_h_gate (per_second).
 * ALGEBRAIC[15] is beta_h in component fast_sodium_current_h_gate (per_second).
 * ALGEBRAIC[37] is i_siCa in component second_inward_current (nanoA).
 * ALGEBRAIC[38] is i_siK in component second_inward_current (nanoA).
 * ALGEBRAIC[40] is i_siNa in component second_inward_current (nanoA).
 * CONSTANTS[28] is P_si in component second_inward_current (nanoA_per_millimolar).
 * STATES[10] is d in component second_inward_current_d_gate (dimensionless).
 * STATES[11] is f in component second_inward_current_f_gate (dimensionless).
 * STATES[12] is f2 in component second_inward_current_f2_gate (dimensionless).
 * ALGEBRAIC[16] is alpha_d in component second_inward_current_d_gate (per_second).
 * ALGEBRAIC[21] is beta_d in component second_inward_current_d_gate (per_second).
 * CONSTANTS[29] is delta_d in component second_inward_current_d_gate (millivolt).
 * ALGEBRAIC[6] is E0_d in component second_inward_current_d_gate (millivolt).
 * ALGEBRAIC[17] is alpha_f in component second_inward_current_f_gate (per_second).
 * ALGEBRAIC[22] is beta_f in component second_inward_current_f_gate (per_second).
 * CONSTANTS[30] is delta_f in component second_inward_current_f_gate (millivolt).
 * ALGEBRAIC[7] is E0_f in component second_inward_current_f_gate (millivolt).
 * CONSTANTS[31] is alpha_f2 in component second_inward_current_f2_gate (per_second).
 * ALGEBRAIC[8] is beta_f2 in component second_inward_current_f2_gate (per_second).
 * CONSTANTS[32] is K_mf2 in component second_inward_current_f2_gate (millimolar).
 * CONSTANTS[33] is radius in component intracellular_sodium_concentration (micrometre).
 * CONSTANTS[34] is length in component intracellular_sodium_concentration (micrometre).
 * CONSTANTS[35] is V_e_ratio in component intracellular_sodium_concentration (dimensionless).
 * CONSTANTS[46] is V_Cell in component intracellular_sodium_concentration (micrometre3).
 * CONSTANTS[47] is Vi in component intracellular_sodium_concentration (micrometre3).
 * CONSTANTS[48] is V_up in component intracellular_calcium_concentration (micrometre3).
 * CONSTANTS[49] is V_rel in component intracellular_calcium_concentration (micrometre3).
 * ALGEBRAIC[39] is i_up in component intracellular_calcium_concentration (nanoA).
 * ALGEBRAIC[41] is i_tr in component intracellular_calcium_concentration (nanoA).
 * ALGEBRAIC[44] is i_rel in component intracellular_calcium_concentration (nanoA).
 * STATES[13] is Ca_up in component intracellular_calcium_concentration (millimolar).
 * STATES[14] is Ca_rel in component intracellular_calcium_concentration (millimolar).
 * CONSTANTS[36] is Ca_up_max in component intracellular_calcium_concentration (millimolar).
 * CONSTANTS[37] is K_mCa in component intracellular_calcium_concentration (millimolar).
 * STATES[15] is p in component intracellular_calcium_concentration (dimensionless).
 * ALGEBRAIC[9] is alpha_p in component intracellular_calcium_concentration (per_second).
 * ALGEBRAIC[18] is beta_p in component intracellular_calcium_concentration (per_second).
 * CONSTANTS[38] is tau_up in component intracellular_calcium_concentration (second).
 * CONSTANTS[39] is tau_rep in component intracellular_calcium_concentration (second).
 * CONSTANTS[40] is tau_rel in component intracellular_calcium_concentration (second).
 * CONSTANTS[41] is rCa in component intracellular_calcium_concentration (dimensionless).
 * CONSTANTS[42] is Ve in component extracellular_potassium_concentration (micrometre3).
 * CONSTANTS[43] is Kb in component extracellular_potassium_concentration (millimolar).
 * ALGEBRAIC[42] is i_mK in component extracellular_potassium_concentration (nanoA).
 * CONSTANTS[44] is pf in component extracellular_potassium_concentration (per_second).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[4] is d/dt y in component hyperpolarising_activated_current_y_gate (dimensionless).
 * RATES[5] is d/dt x in component time_dependent_potassium_current_x_gate (dimensionless).
 * RATES[7] is d/dt s in component transient_outward_current_s_gate (dimensionless).
 * RATES[8] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
 * RATES[9] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
 * RATES[10] is d/dt d in component second_inward_current_d_gate (dimensionless).
 * RATES[11] is d/dt f in component second_inward_current_f_gate (dimensionless).
 * RATES[12] is d/dt f2 in component second_inward_current_f2_gate (dimensionless).
 * RATES[3] is d/dt Nai in component intracellular_sodium_concentration (millimolar).
 * RATES[15] is d/dt p in component intracellular_calcium_concentration (dimensionless).
 * RATES[13] is d/dt Ca_up in component intracellular_calcium_concentration (millimolar).
 * RATES[14] is d/dt Ca_rel in component intracellular_calcium_concentration (millimolar).
 * RATES[6] is d/dt Cai in component intracellular_calcium_concentration (millimolar).
 * RATES[1] is d/dt Kc in component extracellular_potassium_concentration (millimolar).
 * RATES[2] is d/dt Ki in component intracellular_potassium_concentration (millimolar).
 */