#include <assert.h>
#include <stdlib.h>
#include "ten_tusscher_2006.h"


GET_CELL_MODEL_DATA(init_cell_model_data) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;

}

//TODO: this should be called only once for the whole mesh, like in the GPU code
SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    sv[0] = -85.23f;   // V;       millivolt
    sv[1] = 0.00621f;  // Xr1;     dimensionless
    sv[2] = 0.4712f;   // Xr2;     dimensionless
    sv[3] = 0.0095f;   // Xs;      dimensionless
    sv[4] = 0.00172f;  // m;       dimensionless
    sv[5] = 0.7444f;   // h;       dimensionless
    sv[6] = 0.7045f;   // j;       dimensionless
    sv[7] = 3.373e-5f; // d;       dimensionless
    sv[8] = 0.7888f;   // f;       dimensionless
    sv[9] = 0.9755f;   // f2;      dimensionless
    sv[10] = 0.9953f;   // fCass;   dimensionless
    sv[11] = 0.999998f; // s;       dimensionless
    sv[12] = 2.42e-8f;  // r;       dimensionless
    sv[13] = 0.000126f; // Ca_i;    millimolar
    sv[14] = 3.64f;     // Ca_SR;   millimolar
    sv[15] = 0.00036f;  // Ca_ss;   millimolar
    sv[16] = 0.9073f;   // R_prime; dimensionless
    sv[17] = 8.604f;    // Na_i;    millimolar
    sv[18] = 136.89f;   // K_i;     millimolar
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    uint32_t sv_id;
	int i;

    #pragma omp parallel for private(sv_id)
    for (i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);

        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt);

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
    sv[10] = rDY[10];
    sv[11] = rDY[11];
    sv[12] = rDY[12];


    sv[13] = dt*rDY[13] + rY[13];
    sv[14] = dt*rDY[14] + rY[14];
    sv[15] = dt*rDY[15] + rY[15];
    sv[16] = dt*rDY[16] + rY[16];
    sv[17] = dt*rDY[17] + rY[17];
    sv[18] = dt*rDY[18] + rY[18];
}


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {

    // State variables
    const real V = sv[0];      // Membrane variable
    const real Xr1 = sv[1];    // Rapid time dependent potassium current Xr1
    const real Xr2 = sv[2];    // Rapid time dependent potassium current Xr2
    const real Xs = sv[3];     // Slow time dependent potassium current Xs
    const real m = sv[4];      // Fast sodium current m
    const real h = sv[5];      // Fast sodium current h gate
    const real j = sv[6];      // Fast sodium current j gate
    const real d = sv[7];      // L type Ca current d gate
    const real f = sv[8];      // var_L_type_Ca_current_f_gate__f
    const real f2 = sv[9];     // var_L_type_Ca_current_f2_gate__f2
    const real fCass = sv[10]; // L_type_Ca_current__fCass
    const real s = sv[11];     // gating s
    const real r = sv[12];     // gating r
    const real Ca_i = sv[13];  // calcium_dynamics__Ca_i
    const real Ca_SR = sv[14];
    const real Ca_ss = sv[15];
    const real R_prime = sv[16];
    const real Na_i = sv[17]; // var_sodium_dynamics__Na_i
    const real K_i = sv[18];  // var_potassium_dynamics__K_i

    // Some constants
    const real R   = 8314.472;
    const real T   = 310.0;
    const real F   = 96485.3415f;
    const real Cm  = 0.185;
    //const real Cm  = 0.00000001;
    const real V_c = 0.016404;

    const real Ko  = 5.4;
    const real Nao = 140.0;
    const real Cao = 2.0;

    const real P_kna = 0.03;
    const real K_mk  = 1.0;
    const real P_NaK = 2.724;
    const real K_mNa = 40.0;
    const real K_pCa = 0.0005;

    // Calcium dynamics
    const real V_rel    = 0.102;
    const real k1_prime = 0.15;
    const real max_sr   = 2.5;
    const real min_sr   = 1.0;
    const real EC       = 1.5;
    const real Vmax_up  = 0.006375;

    // NCX consts
    const real alpha  = 2.5;
    const real gamma  = 0.35;
    const real K_sat  = 0.1;
    const real Km_Ca  = 1.38;
    const real Km_Nai = 87.5;
    const real K_NaCa = 1000.0;

    const real g_to  = 0.294;
    const real g_Kr  = 0.153;
    const real g_Ks  = 0.098;
    const real g_CaL = 3.98e-05;
    const real g_Na  = 14.838;
    const real g_pK  = 0.0146;
    const real g_bca = 0.000592;
    const real g_pCa = 0.1238;
    const real g_K1  = 5.405;
    const real g_bna = 0.00029;

    // Calculations
    real EK  = ((R * T) / F) * logf(Ko / K_i);
    real EKs = ((R * T) / F) * logf((Ko + (P_kna * Nao)) / (K_i + (P_kna * Na_i)));
    real ENa = ((R * T) / F) * logf(Nao / Na_i);
    real ECa = ((0.5f * R * T) / F) * logf(Cao / Ca_i);

    real beta_K1 = ((3.0f * expf(0.0002f * ((V - EK) + 100.0f))) + expf(0.1f * ((V - EK) - 10.0f))) / (1.0f + expf((-0.5f) * (V - EK)));
    real alpha_K1 = 0.1f / (1.0f + expf(0.06f * ((V - EK) - 200.0f)));
    real xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);

    real IK1 = g_K1 * xK1_inf * (V - EK);
    real Ito = g_to * r * s * (V - EK);
    real IKr = g_Kr * Xr1 * Xr2 * (V - EK) * sqrtf(Ko / 5.4f);
    real IKs = g_Ks * powf(Xs, 2.0f) * (V - EKs);
    real IpK = (g_pK * (V - EK)) / (1.0f + expf((25.0f - V) / 5.98f));

    real ICaL = (V < 15.0f-1.0e-5f || V > 15.0f+1.0e-5f) ? ((((g_CaL * d * f * f2 * fCass * 4.0f * (V - 15.0f) * powf(F, 2.0f)) / (R * T)) * ((0.25f * Ca_ss * expf((2.0f * (V - 15.0f) * F) / (R * T))) - Cao)) / (expf((2.0f * (V - 15.0f) * F) / (R * T)) - 1.0f)) : g_CaL * d * f * f2 * fCass * 2.0f * F * (0.25f * Ca_ss  - Cao);
    real IbCa = g_bca * (V - ECa);
    real IpCa = (g_pCa * Ca_i) / (Ca_i + K_pCa);

    real INaK = ((((P_NaK * Ko) / (Ko + K_mk)) * Na_i) / (Na_i + K_mNa)) / (1.0f + (0.1245f * expf(((-0.1f) * V * F) / (R * T))) + (0.0353f * expf(((-V) * F) / (R * T))));
    real INa  = g_Na * powf(m, 3.0f) * h * j * (V - ENa);
    real IbNa = g_bna * (V - ENa);
    real INaCa = (K_NaCa * ((expf((gamma * V * F) / (R * T)) * powf(Na_i, 3.0f) * Cao) - (expf(((gamma - 1.0f) * V * F) / (R * T)) * powf(Nao, 3) * Ca_i * alpha))) / ((powf(Km_Nai, 3.0f) + powf(Nao, 3.0f)) * (Km_Ca + Cao) * (1.0f + (K_sat * expf(((gamma - 1.0f) * V * F) / (R * T)))));

    // Stimulus
    real var_membrane__i_Stim = stim_current;

    real xr1_inf   = 1.0f / (1.0f + expf(((-26.0f) - V) / 7.0f));
    real alpha_xr1 = 450.0f / (1.0f + expf(((-45.0f) - V) / 10.0f));
    real beta_xr1  = 6.0f / (1.0f + expf((V + 30.0f) / 11.5f));
    real tau_xr1   = 1.0f * alpha_xr1 * beta_xr1;

    real xr2_inf   = 1.0f / (1.0f + expf((V + 88.0f) / 24.0f));
    real alpha_xr2 = 3.0f / (1.0f + expf(((-60.0f) - V) / 20.0f));
    real beta_xr2  = 1.12f / (1.0f + expf((V - 60.0f) / 20.0f));
    real tau_xr2   = 1.0f * alpha_xr2 * beta_xr2;

    real xs_inf   = 1.0f / (1.0f + expf(((-5.0f) - V) / 14.0f));
    real alpha_xs = 1400.0f / sqrtf(1.0f + expf((5.0f - V) / 6.0f));
    real beta_xs  = 1.0f / (1.0f + expf((V - 35.0f) / 15.0f));
    real tau_xs   = (1.0f * alpha_xs * beta_xs) + 80.0f;

    real m_inf   = 1.0f / powf(1.0f + expf(((-56.86f) - V) / 9.03f), 2.0f);
    real alpha_m = 1.0f / (1.0f + expf(((-60.0f) - V) / 5.0f));
    real beta_m  = (0.1f / (1.0f + expf((V + 35.0f) / 5.0f))) + (0.1f / (1.0f + expf((V - 50.0f) / 200.0f)));
    real tau_m   = 1.0f * alpha_m * beta_m;

    real h_inf   = 1.0f / powf(1.0f + expf((V + 71.55f) / 7.43f), 2.0f);
    real alpha_h = (V < (-40.0f)) ? (0.057f * expf((-(V + 80.0f)) / 6.8f)) : 0.0f;
    real beta_h  = (V < (-40.0f)) ? ((2.7f * expf(0.079f * V)) + (310000.0f * expf(0.3485f * V))) : (0.77f / (0.13f * (1.0f + expf((V + 10.66f) / (-11.1f)))));
    real tau_h   = 1.0f / (alpha_h + beta_h);

    real j_inf   = 1.0f / powf(1.0f + expf((V + 71.55f) / 7.43f), 2.0f);
    real alpha_j = (V < (-40.0f)) ? ((((((-25428.0f) * expf(0.2444f * V)) - (6.948e-06f * expf((-0.04391f) * V))) * (V + 37.78f)) / 1.0f) / (1.0f + expf(0.311f * (V + 79.23f)))) : 0.0f;
    real beta_j  = (V < (-40.0f)) ? ((0.02424f * expf((-0.01052f) * V)) / (1.0f + expf((-0.1378f) * (V + 40.14f)))) : ((0.6f * expf(0.057f * V)) / (1.0f + expf((-0.1f) * (V + 32.0f))));
    real tau_j   = 1.0f / (alpha_j + beta_j);

    real d_inf = 1.0f / (1.0f + expf(((-8.0f) - V) / 7.5f));
    real alpha_d = (1.4f / (1.0f + expf(((-35.0f) - V) / 13.0f))) + 0.25f;
    real beta_d  = 1.4f / (1.0f + expf((V + 5.0f) / 5.0f));
    real gamma_d = 1.0f / (1.0f + expf((50.0f - V) / 20.0f));
    real tau_d   = (1.0f * alpha_d * beta_d) + gamma_d;

    real f_inf = 1.0f / (1.0f + expf((V + 20.0f) / 7.0f));
    real tau_f = (1102.5f * expf((-powf(V + 27.0f, 2.0f)) / 225.0f)) + (200.0f / (1.0f + expf((13.0f - V) / 10.0f))) + (180.0f / (1.0f + expf((V + 30.0f) / 10.0f))) + 20.0f;

    real f2_inf = (0.67f / (1.0f + expf((V + 35.0f) / 7.0f))) + 0.33f;
    real tau_f2 = (562.0f * expf((-powf(V + 27.0f, 2.0f)) / 240.0f)) + (31.0f / (1.0f + expf((25.0f - V) / 10.0f))) + (80.0f / (1.0f + expf((V + 30.0f) / 10.0f)));

    real fCass_inf = (0.6f / (1.0f + powf(Ca_ss / 0.05f, 2.0f))) + 0.4f;
    real tau_fCass = (80.0f / (1.0f + powf(Ca_ss / 0.05f, 2.0f))) + 2.0f;

    real s_inf = 1.0f / (1.0f + expf((V + 20.0f) / 5.0f));
    real tau_s = (85.0f * expf((-powf(V + 45.0f, 2.0f)) / 320.0f)) + (5.0f / (1.0f + expf((V - 20.0f) / 5.0f))) + 3.0f;

    real r_inf = 1.0f / (1.0f + expf((20.0f - V) / 6.0f));
    real tau_r = (9.5f * expf((-powf(V + 40.0f, 2.0f)) / 1800.0f)) + 0.8f;

    real kcasr = max_sr - ((max_sr - min_sr) / (1.0f + powf(EC / Ca_SR, 2.0f)));
    real k1 = k1_prime / kcasr;
    const real k3 = 0.06;
    real var_calcium_dynamics__O = (k1 * powf(Ca_ss, 2.0f) * R_prime) / (k3 + (k1 * powf(Ca_ss, 2.0f)));
    real Irel = V_rel * var_calcium_dynamics__O * (Ca_SR - Ca_ss);

    const real var_calcium_dynamics__K_up = 0.00025;
    real var_calcium_dynamics__i_up = Vmax_up / (1.0f + (powf(var_calcium_dynamics__K_up, 2.0f) / powf(Ca_i, 2.0f)));
    real var_calcium_dynamics__V_leak = 0.00036f;
    real var_calcium_dynamics__i_leak = var_calcium_dynamics__V_leak * (Ca_SR - Ca_i);
    const real var_calcium_dynamics__V_xfer = 0.0038f;
    real var_calcium_dynamics__i_xfer = var_calcium_dynamics__V_xfer * (Ca_ss - Ca_i);
    const real var_calcium_dynamics__k2_prime = 0.045f;
    real var_calcium_dynamics__k2 = var_calcium_dynamics__k2_prime * kcasr;
    const real var_calcium_dynamics__k4 = 0.005f;
    const real var_calcium_dynamics__Buf_c = 0.2f;
    const real var_calcium_dynamics__K_buf_c = 0.001f;
    real Ca_i_bufc = 1.0f / (1.0f + ((var_calcium_dynamics__Buf_c * var_calcium_dynamics__K_buf_c) / powf(Ca_i + var_calcium_dynamics__K_buf_c, 2)));
    const real var_calcium_dynamics__K_buf_sr = 0.3f;
    const real var_calcium_dynamics__Buf_sr = 10.0f;
    real var_calcium_dynamics__Ca_sr_bufsr = 1.0f / (1.0f + ((var_calcium_dynamics__Buf_sr * var_calcium_dynamics__K_buf_sr) / powf(Ca_SR + var_calcium_dynamics__K_buf_sr, 2)));
    const real var_calcium_dynamics__Buf_ss = 0.4f;
    const real var_calcium_dynamics__K_buf_ss = 0.00025f;
    real Ca_ss_bufss = 1.0f / (1.0f + ((var_calcium_dynamics__Buf_ss * var_calcium_dynamics__K_buf_ss) / powf(Ca_ss + var_calcium_dynamics__K_buf_ss, 2.0f)));
    const real var_calcium_dynamics__V_sr = 0.001094f;
    const real var_calcium_dynamics__V_ss = 5.468e-05f;
    real var_calcium_dynamics__V_c = V_c;
    real var_calcium_dynamics__F = F;
    real var_calcium_dynamics__Cm = Cm;
    real var_calcium_dynamics__ICaL = ICaL;
    real var_calcium_dynamics__INaCa = INaCa;
    real var_calcium_dynamics__IpCa = IpCa;
    real var_calcium_dynamics__IbCa = IbCa;

    real d_dt_V = -(IK1 + Ito + IKr + IKs + ICaL + INaK + INa + IbNa + INaCa + IbCa + IpK + IpCa + var_membrane__i_Stim);

    real d_dt_R_prime = ((-var_calcium_dynamics__k2) * Ca_ss * R_prime) + (var_calcium_dynamics__k4 * (1.0f - R_prime));
    real d_dt_Ca_i = Ca_i_bufc * (((((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_c) + var_calcium_dynamics__i_xfer) - ((((var_calcium_dynamics__IbCa + var_calcium_dynamics__IpCa) - (2.0f * var_calcium_dynamics__INaCa)) * var_calcium_dynamics__Cm) / (2.0f * var_calcium_dynamics__V_c * var_calcium_dynamics__F)));
    real d_dt_Ca_SR = var_calcium_dynamics__Ca_sr_bufsr * (var_calcium_dynamics__i_up - (Irel + var_calcium_dynamics__i_leak));
    real d_dt_Ca_ss = Ca_ss_bufss * (((((-var_calcium_dynamics__ICaL) * var_calcium_dynamics__Cm) / (2.0f * var_calcium_dynamics__V_ss * var_calcium_dynamics__F)) + ((Irel * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_ss)) - ((var_calcium_dynamics__i_xfer * var_calcium_dynamics__V_c) / var_calcium_dynamics__V_ss));
    real d_dt_Na_i = ((-(INa + IbNa + (3.0f * INaK) + (3.0f * INaCa))) / (V_c * F)) * Cm;
    real d_dt_K_i = ((-((IK1 + Ito + IKr + IKs + IpK + var_membrane__i_Stim) - (2.0f * INaK))) / (V_c * F)) * Cm;

    rDY_[ 0] = d_dt_V;
    rDY_[13] = d_dt_Ca_i;
    rDY_[14] = d_dt_Ca_SR;
    rDY_[15] = d_dt_Ca_ss;
    rDY_[16] = d_dt_R_prime;
    rDY_[17] = d_dt_Na_i;
    rDY_[18] = d_dt_K_i;


    // Rush Larsen
    rDY_[ 1] = xr1_inf + (Xr1-xr1_inf)*expf(-dt/tau_xr1) ;
    rDY_[ 2] = xr2_inf + (Xr2-xr2_inf)*expf(-dt/tau_xr2);
    rDY_[ 3] = xs_inf + (Xs-xs_inf)*expf(-dt/tau_xs);
    rDY_[ 4] = m_inf + (m-m_inf)*expf(-dt/tau_m);
    rDY_[ 5] = h_inf + (h-h_inf)*expf(-dt/tau_h);
    rDY_[ 6] = j_inf + (j-j_inf)*expf(-dt/tau_j);
    rDY_[ 7] = d_inf + (d-d_inf)*expf(-dt/tau_d);
    rDY_[ 8] = f_inf + (f-f_inf)*expf(-dt/tau_f);
    rDY_[ 9] = f2_inf + (f2-f2_inf)*expf(-dt/tau_f2);
    rDY_[10] = fCass_inf + (fCass-fCass_inf)*expf(-dt/tau_fCass);
    rDY_[11] = s_inf + (s-s_inf)*expf(-dt/tau_s);
    rDY_[12] = r_inf + (r-r_inf)*expf(-dt/tau_r);


}