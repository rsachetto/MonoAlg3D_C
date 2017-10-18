#include "model_common.h"
#include <assert.h>
#include <stdlib.h>
#include <unitypes.h>

void RHS_cpu(const Real *sv, Real *rDY_, Real stim_current, Real dt);
void solve_model_ode_cpu(Real dt, Real *sv, Real stim_current, int neq);

        void init_cell_model_data(struct cell_model_data* cell_model, bool get_initial_v, bool get_neq) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = -85.23f;
    if(get_neq)
        cell_model->number_of_ode_equations = 19;

}

void set_model_initial_conditions_cpu(Real *sv) {

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

void solve_model_odes_cpu(Real dt, Real *sv, Real *stim_currents, uint32_t *cells_to_solve,
                          uint32_t num_cells_to_solve,int num_steps, int neq, void *extra_data) {

    uint32_t sv_id;

    #pragma omp parallel for private(sv_id)
    for (int i = 0; i < num_cells_to_solve; i++) {
        sv_id = cells_to_solve[i];

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * neq), stim_currents[i], neq);

        }
    }
}

void solve_model_ode_cpu(Real dt, Real *sv, Real stim_current, int neq)  {

    assert(sv);

    Real rY[neq], rDY[neq];

    for(int i = 0; i < neq; i++)
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


void RHS_cpu(const Real *sv, Real *rDY_, Real stim_current, Real dt) {

    // State variables
    const Real V = sv[0];      // Membrane variable
    const Real Xr1 = sv[1];    // Rapid time dependent potassium current Xr1
    const Real Xr2 = sv[2];    // Rapid time dependent potassium current Xr2
    const Real Xs = sv[3];     // Slow time dependent potassium current Xs
    const Real m = sv[4];      // Fast sodium current m
    const Real h = sv[5];      // Fast sodium current h gate
    const Real j = sv[6];      // Fast sodium current j gate
    const Real d = sv[7];      // L type Ca current d gate
    const Real f = sv[8];      // var_L_type_Ca_current_f_gate__f
    const Real f2 = sv[9];     // var_L_type_Ca_current_f2_gate__f2
    const Real fCass = sv[10]; // L_type_Ca_current__fCass
    const Real s = sv[11];     // gating s
    const Real r = sv[12];     // gating r
    const Real Ca_i = sv[13];  // calcium_dynamics__Ca_i
    const Real Ca_SR = sv[14];
    const Real Ca_ss = sv[15];
    const Real R_prime = sv[16];
    const Real Na_i = sv[17]; // var_sodium_dynamics__Na_i
    const Real K_i = sv[18];  // var_potassium_dynamics__K_i

    // Some constants
    const Real R   = 8314.472;
    const Real T   = 310.0;
    const Real F   = 96485.3415f;
    const Real Cm  = 0.185;
    //const Real Cm  = 0.00000001;
    const Real V_c = 0.016404;

    const Real Ko  = 5.4;
    const Real Nao = 140.0;
    const Real Cao = 2.0;

    const Real P_kna = 0.03;
    const Real K_mk  = 1.0;
    const Real P_NaK = 2.724;
    const Real K_mNa = 40.0;
    const Real K_pCa = 0.0005;

    // Calcium dynamics
    const Real V_rel    = 0.102;
    const Real k1_prime = 0.15;
    const Real max_sr   = 2.5;
    const Real min_sr   = 1.0;
    const Real EC       = 1.5;
    const Real Vmax_up  = 0.006375;

    // NCX consts
    const Real alpha  = 2.5;
    const Real gamma  = 0.35;
    const Real K_sat  = 0.1;
    const Real Km_Ca  = 1.38;
    const Real Km_Nai = 87.5;
    const Real K_NaCa = 1000.0;

    const Real g_to  = 0.294;
    const Real g_Kr  = 0.153;
    const Real g_Ks  = 0.098;
    const Real g_CaL = 3.98e-05;
    const Real g_Na  = 14.838;
    const Real g_pK  = 0.0146;
    const Real g_bca = 0.000592;
    const Real g_pCa = 0.1238;
    const Real g_K1  = 5.405;
    const Real g_bna = 0.00029;

    // Calculations
    Real EK  = ((R * T) / F) * logf(Ko / K_i);
    Real EKs = ((R * T) / F) * logf((Ko + (P_kna * Nao)) / (K_i + (P_kna * Na_i)));
    Real ENa = ((R * T) / F) * logf(Nao / Na_i);
    Real ECa = ((0.5f * R * T) / F) * logf(Cao / Ca_i);

    Real beta_K1 = ((3.0f * expf(0.0002f * ((V - EK) + 100.0f))) + expf(0.1f * ((V - EK) - 10.0f))) / (1.0f + expf((-0.5f) * (V - EK)));
    Real alpha_K1 = 0.1f / (1.0f + expf(0.06f * ((V - EK) - 200.0f)));
    Real xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);

    Real IK1 = g_K1 * xK1_inf * (V - EK);
    Real Ito = g_to * r * s * (V - EK);
    Real IKr = g_Kr * Xr1 * Xr2 * (V - EK) * sqrtf(Ko / 5.4f);
    Real IKs = g_Ks * powf(Xs, 2.0f) * (V - EKs);
    Real IpK = (g_pK * (V - EK)) / (1.0f + expf((25.0f - V) / 5.98f));

    Real ICaL = (V < 15.0f-1.0e-5f || V > 15.0f+1.0e-5f) ? ((((g_CaL * d * f * f2 * fCass * 4.0f * (V - 15.0f) * powf(F, 2.0f)) / (R * T)) * ((0.25f * Ca_ss * expf((2.0f * (V - 15.0f) * F) / (R * T))) - Cao)) / (expf((2.0f * (V - 15.0f) * F) / (R * T)) - 1.0f)) : g_CaL * d * f * f2 * fCass * 2.0f * F * (0.25f * Ca_ss  - Cao);
    Real IbCa = g_bca * (V - ECa);
    Real IpCa = (g_pCa * Ca_i) / (Ca_i + K_pCa);

    Real INaK = ((((P_NaK * Ko) / (Ko + K_mk)) * Na_i) / (Na_i + K_mNa)) / (1.0f + (0.1245f * expf(((-0.1f) * V * F) / (R * T))) + (0.0353f * expf(((-V) * F) / (R * T))));
    Real INa  = g_Na * powf(m, 3.0f) * h * j * (V - ENa);
    Real IbNa = g_bna * (V - ENa);
    Real INaCa = (K_NaCa * ((expf((gamma * V * F) / (R * T)) * powf(Na_i, 3.0f) * Cao) - (expf(((gamma - 1.0f) * V * F) / (R * T)) * powf(Nao, 3) * Ca_i * alpha))) / ((powf(Km_Nai, 3.0f) + powf(Nao, 3.0f)) * (Km_Ca + Cao) * (1.0f + (K_sat * expf(((gamma - 1.0f) * V * F) / (R * T)))));

    // Stimulus
    Real var_membrane__i_Stim = stim_current;

    Real xr1_inf   = 1.0f / (1.0f + expf(((-26.0f) - V) / 7.0f));
    Real alpha_xr1 = 450.0f / (1.0f + expf(((-45.0f) - V) / 10.0f));
    Real beta_xr1  = 6.0f / (1.0f + expf((V + 30.0f) / 11.5f));
    Real tau_xr1   = 1.0f * alpha_xr1 * beta_xr1;

    Real xr2_inf   = 1.0f / (1.0f + expf((V + 88.0f) / 24.0f));
    Real alpha_xr2 = 3.0f / (1.0f + expf(((-60.0f) - V) / 20.0f));
    Real beta_xr2  = 1.12f / (1.0f + expf((V - 60.0f) / 20.0f));
    Real tau_xr2   = 1.0f * alpha_xr2 * beta_xr2;

    Real xs_inf   = 1.0f / (1.0f + expf(((-5.0f) - V) / 14.0f));
    Real alpha_xs = 1400.0f / sqrtf(1.0f + expf((5.0f - V) / 6.0f));
    Real beta_xs  = 1.0f / (1.0f + expf((V - 35.0f) / 15.0f));
    Real tau_xs   = (1.0f * alpha_xs * beta_xs) + 80.0f;

    Real m_inf   = 1.0f / powf(1.0f + expf(((-56.86f) - V) / 9.03f), 2.0f);
    Real alpha_m = 1.0f / (1.0f + expf(((-60.0f) - V) / 5.0f));
    Real beta_m  = (0.1f / (1.0f + expf((V + 35.0f) / 5.0f))) + (0.1f / (1.0f + expf((V - 50.0f) / 200.0f)));
    Real tau_m   = 1.0f * alpha_m * beta_m;

    Real h_inf   = 1.0f / powf(1.0f + expf((V + 71.55f) / 7.43f), 2.0f);
    Real alpha_h = (V < (-40.0f)) ? (0.057f * expf((-(V + 80.0f)) / 6.8f)) : 0.0f;
    Real beta_h  = (V < (-40.0f)) ? ((2.7f * expf(0.079f * V)) + (310000.0f * expf(0.3485f * V))) : (0.77f / (0.13f * (1.0f + expf((V + 10.66f) / (-11.1f)))));
    Real tau_h   = 1.0f / (alpha_h + beta_h);

    Real j_inf   = 1.0f / powf(1.0f + expf((V + 71.55f) / 7.43f), 2.0f);
    Real alpha_j = (V < (-40.0f)) ? ((((((-25428.0f) * expf(0.2444f * V)) - (6.948e-06f * expf((-0.04391f) * V))) * (V + 37.78f)) / 1.0f) / (1.0f + expf(0.311f * (V + 79.23f)))) : 0.0f;
    Real beta_j  = (V < (-40.0f)) ? ((0.02424f * expf((-0.01052f) * V)) / (1.0f + expf((-0.1378f) * (V + 40.14f)))) : ((0.6f * expf(0.057f * V)) / (1.0f + expf((-0.1f) * (V + 32.0f))));
    Real tau_j   = 1.0f / (alpha_j + beta_j);

    Real d_inf = 1.0f / (1.0f + expf(((-8.0f) - V) / 7.5f));
    Real alpha_d = (1.4f / (1.0f + expf(((-35.0f) - V) / 13.0f))) + 0.25f;
    Real beta_d  = 1.4f / (1.0f + expf((V + 5.0f) / 5.0f));
    Real gamma_d = 1.0f / (1.0f + expf((50.0f - V) / 20.0f));
    Real tau_d   = (1.0f * alpha_d * beta_d) + gamma_d;

    Real f_inf = 1.0f / (1.0f + expf((V + 20.0f) / 7.0f));
    Real tau_f = (1102.5f * expf((-powf(V + 27.0f, 2.0f)) / 225.0f)) + (200.0f / (1.0f + expf((13.0f - V) / 10.0f))) + (180.0f / (1.0f + expf((V + 30.0f) / 10.0f))) + 20.0f;

    Real f2_inf = (0.67f / (1.0f + expf((V + 35.0f) / 7.0f))) + 0.33f;
    Real tau_f2 = (562.0f * expf((-powf(V + 27.0f, 2.0f)) / 240.0f)) + (31.0f / (1.0f + expf((25.0f - V) / 10.0f))) + (80.0f / (1.0f + expf((V + 30.0f) / 10.0f)));

    Real fCass_inf = (0.6f / (1.0f + powf(Ca_ss / 0.05f, 2.0f))) + 0.4f;
    Real tau_fCass = (80.0f / (1.0f + powf(Ca_ss / 0.05f, 2.0f))) + 2.0f;

    Real s_inf = 1.0f / (1.0f + expf((V + 20.0f) / 5.0f));
    Real tau_s = (85.0f * expf((-powf(V + 45.0f, 2.0f)) / 320.0f)) + (5.0f / (1.0f + expf((V - 20.0f) / 5.0f))) + 3.0f;

    Real r_inf = 1.0f / (1.0f + expf((20.0f - V) / 6.0f));
    Real tau_r = (9.5f * expf((-powf(V + 40.0f, 2.0f)) / 1800.0f)) + 0.8f;

    Real kcasr = max_sr - ((max_sr - min_sr) / (1.0f + powf(EC / Ca_SR, 2.0f)));
    Real k1 = k1_prime / kcasr;
    const Real k3 = 0.06;
    Real var_calcium_dynamics__O = (k1 * powf(Ca_ss, 2.0f) * R_prime) / (k3 + (k1 * powf(Ca_ss, 2.0f)));
    Real Irel = V_rel * var_calcium_dynamics__O * (Ca_SR - Ca_ss);

    const Real var_calcium_dynamics__K_up = 0.00025;
    Real var_calcium_dynamics__i_up = Vmax_up / (1.0f + (powf(var_calcium_dynamics__K_up, 2.0f) / powf(Ca_i, 2.0f)));
    Real var_calcium_dynamics__V_leak = 0.00036f;
    Real var_calcium_dynamics__i_leak = var_calcium_dynamics__V_leak * (Ca_SR - Ca_i);
    const Real var_calcium_dynamics__V_xfer = 0.0038f;
    Real var_calcium_dynamics__i_xfer = var_calcium_dynamics__V_xfer * (Ca_ss - Ca_i);
    const Real var_calcium_dynamics__k2_prime = 0.045f;
    Real var_calcium_dynamics__k2 = var_calcium_dynamics__k2_prime * kcasr;
    const Real var_calcium_dynamics__k4 = 0.005f;
    const Real var_calcium_dynamics__Buf_c = 0.2f;
    const Real var_calcium_dynamics__K_buf_c = 0.001f;
    Real Ca_i_bufc = 1.0f / (1.0f + ((var_calcium_dynamics__Buf_c * var_calcium_dynamics__K_buf_c) / powf(Ca_i + var_calcium_dynamics__K_buf_c, 2)));
    const Real var_calcium_dynamics__K_buf_sr = 0.3f;
    const Real var_calcium_dynamics__Buf_sr = 10.0f;
    Real var_calcium_dynamics__Ca_sr_bufsr = 1.0f / (1.0f + ((var_calcium_dynamics__Buf_sr * var_calcium_dynamics__K_buf_sr) / powf(Ca_SR + var_calcium_dynamics__K_buf_sr, 2)));
    const Real var_calcium_dynamics__Buf_ss = 0.4f;
    const Real var_calcium_dynamics__K_buf_ss = 0.00025f;
    Real Ca_ss_bufss = 1.0f / (1.0f + ((var_calcium_dynamics__Buf_ss * var_calcium_dynamics__K_buf_ss) / powf(Ca_ss + var_calcium_dynamics__K_buf_ss, 2.0f)));
    const Real var_calcium_dynamics__V_sr = 0.001094f;
    const Real var_calcium_dynamics__V_ss = 5.468e-05f;
    Real var_calcium_dynamics__V_c = V_c;
    Real var_calcium_dynamics__F = F;
    Real var_calcium_dynamics__Cm = Cm;
    Real var_calcium_dynamics__ICaL = ICaL;
    Real var_calcium_dynamics__INaCa = INaCa;
    Real var_calcium_dynamics__IpCa = IpCa;
    Real var_calcium_dynamics__IbCa = IbCa;

    Real d_dt_V = -(IK1 + Ito + IKr + IKs + ICaL + INaK + INa + IbNa + INaCa + IbCa + IpK + IpCa + var_membrane__i_Stim);

    Real d_dt_R_prime = ((-var_calcium_dynamics__k2) * Ca_ss * R_prime) + (var_calcium_dynamics__k4 * (1.0f - R_prime));
    Real d_dt_Ca_i = Ca_i_bufc * (((((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_c) + var_calcium_dynamics__i_xfer) - ((((var_calcium_dynamics__IbCa + var_calcium_dynamics__IpCa) - (2.0f * var_calcium_dynamics__INaCa)) * var_calcium_dynamics__Cm) / (2.0f * var_calcium_dynamics__V_c * var_calcium_dynamics__F)));
    Real d_dt_Ca_SR = var_calcium_dynamics__Ca_sr_bufsr * (var_calcium_dynamics__i_up - (Irel + var_calcium_dynamics__i_leak));
    Real d_dt_Ca_ss = Ca_ss_bufss * (((((-var_calcium_dynamics__ICaL) * var_calcium_dynamics__Cm) / (2.0f * var_calcium_dynamics__V_ss * var_calcium_dynamics__F)) + ((Irel * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_ss)) - ((var_calcium_dynamics__i_xfer * var_calcium_dynamics__V_c) / var_calcium_dynamics__V_ss));
    Real d_dt_Na_i = ((-(INa + IbNa + (3.0f * INaK) + (3.0f * INaCa))) / (V_c * F)) * Cm;
    Real d_dt_K_i = ((-((IK1 + Ito + IKr + IKs + IpK + var_membrane__i_Stim) - (2.0f * INaK))) / (V_c * F)) * Cm;

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