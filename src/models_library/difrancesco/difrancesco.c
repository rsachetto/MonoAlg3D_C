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

    log_info("Using DiFrancesco & Noble 1985 CPU model\n");

    uint32_t num_volumes = solver->original_num_cells;

    solver->sv = MALLOC_ARRAY_OF_TYPE(real, NEQ*num_volumes);

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

    // Pacing BCL = 300 ms
    /*
    sv[0] = -89.523;    // V millivolt
    sv[1] = 4.25602;      // Kc millimolar
    sv[2] = 139.882;    // Ki millimolar
    sv[3] = 8.09539;      // Nai millimolar
    sv[4] = 0.0826828;    // y dimensionless
    sv[5] = 0.194361;   // x dimensionless
    sv[6] = 2.57238e-05;   // Cai millimolar
    sv[7] = 0.420488;      // s dimensionless
    sv[8] = 0.00255133;   // m dimensionless
    sv[9] = 0.986752;    // h dimensionless
    sv[10] = 5.50534e-08; // d dimensionless
    sv[11] = 1.0;     // f dimensionless
    sv[12] = 0.742853;     // f2 dimensionless
    sv[13] = 2.09383;     // Ca_up millimolar
    sv[14] = 0.205441;     // Ca_rel millimolar
    sv[15] = 0.987422;     // p dimensionless
    */

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


    //printf("\nInside\n");
    //for (uint32_t i = 0; i < 10; i++)
    //    printf("%g\n",sv[i]);

/*
    //State variables
    const real V_old_ = sv[0];
    const real kc_old_ = sv[1];
    const real ki_old_ = sv[2];
    const real nai_old_ = sv[3];
    const real y_old_ = sv[4];
    const real x_old_ = sv[5];
    const real cai_old_ = sv[6];
    const real s_old_ = sv[7];
    const real m_old_ = sv[8];
    const real h_old_ = sv[9];
    const real d_old_ = sv[10];
    const real f_old_ = sv[11];
    const real f2_old_ = sv[12];
    const real ca_up_old_ = sv[13];
    const real ca_rel_old_ = sv[14];
    const real p_old_ = sv[15];

    //Parameters
    const real R = 8314.472;
    const real T = 310;
    const real F = 96485.3415;
    const real C = 0.075;
    const real i_pulse = 0;
    const real g_f_Na = 3;
    const real g_f_K = 3;
    const real Km_f = 45;
    const real Nao = 140;
    const real delta_y = 1e-5;
    const real i_K_max = 180;
    const real g_K1 = 920;
    const real Km_K1 = 210;
    const real Km_to = 10;
    const real Km_Ca = 0.0005;
    const real g_to = 0.28;
    const real g_Nab = 0.18;
    const real g_Cab = 0.02;
    const real Cao = 2;
    const real I_p = 125;
    const real K_mK = 1;
    const real K_mNa = 40;
    const real n_NaCa = 3;
    const real K_NaCa = 0.02;
    const real d_NaCa = 0.001;
    const real gamma = 0.5;
    const real g_Na = 750;
    const real delta_m = 1e-5;
    const real P_si = 15;
    const real delta_d = 0.0001;
    const real delta_f = 0.0001;
    const real alpha_f2 = 5;
    const real K_mf2 = 0.001;
    const real radius = 0.05;
    const real length = 2;
    const real V_e_ratio = 0.1;
    const real Ca_up_max  = 5;
    const real K_mCa = 0.001;
    const real tau_up = 0.025;
    const real tau_rep = 2;
    const real tau_rel = 0.05;
    const real rCa = 2;
    const real Ve = 0.00157;
    const real Kb = 4;
    const real pf = 0.7;
    const real RTONF = ( R*T)/F;
    const real V_cell =  3.14159*pow(radius, 2.00000)*length;
    const real Vi =  V_cell*(1.00000 - V_e_ratio);
    const real V_up =  Vi*0.0500000;
    const real V_rel =  Vi*0.0200000;

    // Algebraics
    const real beta_f2 = ( cai_old_*alpha_f2)/K_mf2;
    const real alpha_x = ( 0.500000*exp( 0.0826000*(V_old_+50.0000)))/(1.00000+exp( 0.0570000*(V_old_+50.0000)));
    const real beta_x = ( 1.30000*exp( - 0.0600000*(V_old_+20.0000)))/(1.00000+exp( - 0.0400000*(V_old_+20.0000)));

    const real alpha_s =  0.0330000*exp(- V_old_/17.0000);
    const real beta_s = 33.0000/(1.00000+exp(- (V_old_+10.0000)/8.00000));

    const real alpha_h =  20.0000*exp( - 0.125000*(V_old_+75.0000));
    const real beta_h = 2000.00/( 320.000*exp( - 0.100000*(V_old_+75.0000))+1.00000);

    const real alpha_p = ( 0.625000*(V_old_+34.0000))/(exp((V_old_+34.0000)/4.00000) - 1.00000);
    const real beta_p = 5.00000/(1.00000+exp(( - 1.00000*(V_old_+34.0000))/4.00000));

    const real alpha_y =  0.0500000*exp( - 0.0670000*((V_old_+52.0000) - 10.0000));
    const real E0_y = (V_old_+52.0000) - 10.0000;
    const real beta_y = (fabs(E0_y)<delta_y ? 2.50000 : ( 1.00000*E0_y)/(1.00000 - exp( - 0.200000*E0_y)));

    const real E0_m = V_old_+41.0000;
    const real alpha_m = (fabs(E0_m)<delta_m ? 2000.00 : ( 200.000*E0_m)/(1.00000 - exp( - 0.100000*E0_m)));
    const real beta_m =  8000.00*exp( - 0.0560000*(V_old_+66.0000));

    const real E0_d = (V_old_+24.0000) - 5.00000;
    const real alpha_d = (fabs(E0_d)<delta_d ? 120.000 : ( 30.0000*E0_d)/(1.00000 - exp(( - 1.00000*E0_d)/4.00000)));
    const real beta_d = (fabs(E0_d)<delta_d ? 120.000 : ( 12.0000*E0_d)/(exp(E0_d/10.0000) - 1.00000));

    const real E0_f = V_old_+34.0000;
    const real alpha_f = (fabs(E0_f)<delta_f ? 25.0000 : ( 6.25000*E0_f)/(exp(E0_f/4.00000) - 1.00000));
    const real beta_f = 50.0000/(1.00000+exp(( - 1.00000*(V_old_+34.0000))/4.00000));

    const real E_Na =  RTONF*log(Nao/nai_old_);
    const real i_Na_b =  g_Nab*(V_old_ - E_Na);
    const real i_p = ( (( I_p*kc_old_)/(K_mK+kc_old_))*nai_old_)/(K_mNa+nai_old_);
    const real i_NaCa = ( K_NaCa*( exp(( gamma*(n_NaCa - 2.00000)*V_old_)/RTONF)*pow(nai_old_, n_NaCa)*Cao -  exp(( (gamma - 1.00000)*(n_NaCa - 2.00000)*V_old_)/RTONF)*pow(Nao, n_NaCa)*cai_old_))/( (1.00000+ d_NaCa*( cai_old_*pow(Nao, n_NaCa)+ Cao*pow(nai_old_, n_NaCa)))*(1.00000+cai_old_/0.00690000));
    const real E_mh =  RTONF*log((Nao+ 0.120000*kc_old_)/(nai_old_+ 0.120000*ki_old_));
    const real i_Na =  g_Na*pow(m_old_, 3.00000)*h_old_*(V_old_ - E_mh);
    const real i_fNa =  (( y_old_*kc_old_)/(kc_old_+Km_f))*g_f_Na*(V_old_ - E_Na);
    const real i_siNa =  (( 0.0100000*P_si*(V_old_ - 50.0000))/( RTONF*(1.00000 - exp(( - 1.00000*(V_old_ - 50.0000))/RTONF))))*( nai_old_*exp(50.0000/RTONF) -  Nao*exp(( - 1.00000*(V_old_ - 50.0000))/RTONF))*d_old_*f_old_*f2_old_;

    const real i_up =  (( 2.00000*1.00000*Vi*F)/( 1.00000*tau_up*Ca_up_max))*cai_old_*(Ca_up_max - ca_up_old_);
    const real i_tr =  (( 2.00000*1.00000*V_rel*F)/( 1.00000*tau_rep))*p_old_*(ca_up_old_ - ca_rel_old_);

    const real I_K = ( i_K_max*(ki_old_ -  kc_old_*exp(- V_old_/RTONF)))/140.000;
    const real i_K =  x_old_*I_K;
    const real E_K =  RTONF*log(kc_old_/ki_old_);
    const real i_K1 = ( (( g_K1*kc_old_)/(kc_old_+Km_K1))*(V_old_ - E_K))/(1.00000+exp(( ((V_old_+10.0000) - E_K)*2.00000)/RTONF));
    const real i_to =  (( (( s_old_*g_to*(0.200000+kc_old_/(Km_to+kc_old_))*cai_old_)/(Km_Ca+cai_old_))*(V_old_+10.0000))/(1.00000 - exp( - 0.200000*(V_old_+10.0000))))*( ki_old_*exp(( 0.500000*V_old_)/RTONF) -  kc_old_*exp(( - 0.500000*V_old_)/RTONF));
    const real i_fK =  (( y_old_*kc_old_)/(kc_old_+Km_f))*g_f_K*(V_old_ - E_K);
    const real i_siK =  (( 0.0100000*P_si*(V_old_ - 50.0000))/( RTONF*(1.00000 - exp(( - 1.00000*(V_old_ - 50.0000))/RTONF))))*( ki_old_*exp(50.0000/RTONF) -  kc_old_*exp(( - 1.00000*(V_old_ - 50.0000))/RTONF))*d_old_*f_old_*f2_old_;
    const real i_mK = (i_K1+i_K+i_fK+i_siK+i_to) -  2.00000*i_p;

    const real i_f = i_fNa+i_fK;
    const real E_Ca =  0.500000*RTONF*log(Cao/cai_old_);
    const real i_Ca_b =  g_Cab*(V_old_ - E_Ca);
    const real i_siCa =  (( 4.00000*P_si*(V_old_ - 50.0000))/( RTONF*(1.00000 - exp(( - 1.00000*(V_old_ - 50.0000)*2.00000)/RTONF))))*( cai_old_*exp(100.000/RTONF) -  Cao*exp(( - 2.00000*(V_old_ - 50.0000))/RTONF))*d_old_*f_old_*f2_old_;
    const real i_si = i_siCa+i_siK+i_siNa;

    const real i_rel = ( (( 2.00000*1.00000*V_rel*F)/( 1.00000*tau_rel))*ca_rel_old_*pow(cai_old_, rCa))/(pow(cai_old_, rCa)+pow(K_mCa, rCa));

    // Rates of change
    rDY_[0] = - (i_f+i_K+i_K1+i_to+i_Na_b+i_Ca_b+i_p+i_NaCa+i_Na+i_si+i_pulse)/C;
    rDY_[1] =  - pf*(kc_old_ - Kb)+( 1.00000*i_mK)/( 1.00000*Ve*F);
    rDY_[2] = ( - 1.00000*i_mK)/( 1.00000*Vi*F);
    rDY_[3] = ( - 1.00000*(i_Na+i_Na_b+i_fNa+i_siNa+ i_p*3.00000+( i_NaCa*n_NaCa)/(n_NaCa - 2.00000)))/( 1.00000*Vi*F);
    rDY_[4] =  alpha_y*(1.00000 - y_old_) -  beta_y*y_old_;
    rDY_[5] =  alpha_x*(1.00000 - x_old_) -  beta_x*x_old_;
    rDY_[6] = ( - 1.00000*((((i_siCa+i_Ca_b) - ( 2.00000*i_NaCa)/(n_NaCa - 2.00000)) - i_rel)+i_up))/( 2.00000*1.00000*Vi*F);
    rDY_[7] =  alpha_s*(1.00000 - s_old_) -  beta_s*s_old_;
    rDY_[8] =  alpha_m*(1.00000 - m_old_) -  beta_m*m_old_;
    rDY_[9] =  alpha_h*(1.00000 - h_old_) -  beta_h*h_old_;
    rDY_[10] =  alpha_d*(1.00000 - d_old_) -  beta_d*d_old_;
    rDY_[11] =  alpha_f*(1.00000 - f_old_) -  beta_f*f_old_;
    rDY_[12] = alpha_f2 -  f2_old_*(alpha_f2+beta_f2);
    rDY_[13] = ( 1.00000*(i_up - i_tr))/( 2.00000*1.00000*V_up*F);
    rDY_[14] = ( 1.00000*(i_tr - i_rel))/( 2.00000*1.00000*V_rel*F);
    rDY_[15] =  alpha_p*(1.00000 - p_old_) -  beta_p*p_old_;
*/

}

