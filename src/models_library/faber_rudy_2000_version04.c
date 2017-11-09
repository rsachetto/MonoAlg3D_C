#include "faber_rudy_2000_version04.h"
#include <stdlib.h>
GET_CELL_MODEL_DATA (init_cell_model_data) {

    if (get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if (get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU (set_model_initial_conditions_cpu) {

    sv[0] = -90.000000f;  // V millivolt
    sv[1] = 0.000800f;    // m dimensionless
    sv[2] = 0.993771f;    // h dimensionless
    sv[3] = 0.995727f;    // j dimensionless
    sv[4] = 0.000003f;    // d dimensionless
    sv[5] = 0.999837f;    // f dimensionless
    sv[6] = 0.000970f;    // b dimensionless
    sv[7] = 0.994305f;    // g dimensionless
    sv[8] = 0.000124f;    // xr dimensionless
    sv[9] = 0.004457f;    // xs1 dimensionless
    sv[10] = 0.004457f;   // xs2 dimensionless
    sv[11] = 0.500000f;   // zdv dimensionless
    sv[12] = 0.500000f;   // ydv dimensionless
    sv[13] = 0.000000f;   // APtrack dimensionless
    sv[14] = 0.000000f;   // APtrack2 dimensionless
    sv[15] = 0.000000f;   // APtrack3 dimensionless
    sv[16] = 0.000000f;   // Cainfluxtrack dimensionless
    sv[17] = 0.000000f;   // OVRLDtrack dimensionless
    sv[18] = 0.000000f;   // OVRLDtrack2 dimensionless
    sv[19] = 0.000000f;   // OVRLDtrack3 dimensionless
    sv[20] = 1.800000f;   // Ca_JSR millimolar
    sv[21] = 1.800000f;   // Ca_NSR millimolar
    sv[22] = 0.000060f;   // Cai millimolar
    sv[23] = 9.000000f;   // Nai millimolar
    sv[24] = 141.200000f; // Ki millimolar
}

SOLVE_MODEL_ODES_CPU (solve_model_odes_cpu) {

    uint32_t sv_id;

#pragma omp parallel for private(sv_id)
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if (cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu (dt, sv + (sv_id * NEQ), stim_currents[i]);
        }
    }
}

void solve_model_ode_cpu (real dt, real *sv, real stim_current) {

    real rY[NEQ], rDY[NEQ];

    for (int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu (rY, rDY, stim_current);

    for (int i = 0; i < NEQ; i++)
        sv[i] = dt * rDY[i] + rY[i];
}

#define IFNUMBER_1(name)                                                                                               \
    if ((fabsf (calc_E0_m) >= delta_m)) {                                                                              \
         (name) = ((3.200000000000000e+02f * calc_E0_m) /                                                              \
                   (1.000000000000000e+00f - expf (((-1.000000000000000e-01f) * calc_E0_m))));                         \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = 3.200000000000000e+03f;                                                                               \
    }
#define IFNUMBER_2(name)                                                                                               \
    if ((V_old_ < (-4.000000000000000e+01f))) {                                                                        \
         (name) = (1.350000000000000e+02f * expf (((8.000000000000000e+01f + V_old_) / (-6.800000000000000e+00f))));   \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = 0.000000000000000e+00f;                                                                               \
    }
#define IFNUMBER_3(name)                                                                                               \
    if ((V_old_ < (-4.000000000000000e+01f))) {                                                                        \
         (name) = ((3.560000000000000e+03f * expf ((7.900000000000000e-02f * V_old_))) +                               \
                   (3.100000000000000e+08f * expf ((3.500000000000000e-01f * V_old_))));                               \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) =                                                                                                       \
            (1.000000000000000e+03f /                                                                                  \
             (1.300000000000000e-01f *                                                                                 \
              (1.000000000000000e+00f + expf (((V_old_ + 1.066000000000000e+01f) / (-1.110000000000000e+01f))))));     \
    }
#define IFNUMBER_4(name)                                                                                               \
    if ((V_old_ < (-4.000000000000000e+01f))) {                                                                        \
         (name) = ((1.000000000000000e+03f *                                                                           \
                    (-((1.271400000000000e+05f * expf ((2.444000000000000e-01f * V_old_))) +                           \
                       (3.474000000000000e-05f * expf (((-4.391000000000000e-02f) * V_old_))))) *                      \
                    (V_old_ + 3.778000000000000e+01f)) /                                                               \
                   (1.000000000000000e+00f + expf ((3.110000000000000e-01f * (V_old_ + 7.923000000000000e+01f)))));    \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = 0.000000000000000e+00f;                                                                               \
    }
#define IFNUMBER_5(name)                                                                                               \
    if ((V_old_ < (-4.000000000000000e+01f))) {                                                                        \
         (name) = ((1.212000000000000e+02f * expf (((-1.052000000000000e-02f) * V_old_))) /                            \
                   (1.000000000000000e+00f + expf (((-1.378000000000000e-01f) * (V_old_ + 4.014000000000000e+01f))))); \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((3.000000000000000e+02f * expf (((-2.535000000000000e-07f) * V_old_))) /                             \
                  (1.000000000000000e+00f + expf (((-1.000000000000000e-01f) * (V_old_ + 3.200000000000000e+01f)))));  \
    }
#define IFNUMBER_6(name)                                                                                               \
    if ((fabsf (calc_E0_d) < 1.000000000000000e-05f)) {                                                                \
         (name) = (1.000000000000000e-03f / (3.500000000000000e-02f * 6.240000000000000e+00f));                        \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((1.000000000000000e-03f * calc_d_infinity *                                                          \
                   (1.000000000000000e+00f - expf (((-calc_E0_d) / 6.240000000000000e+00f)))) /                        \
                  (3.500000000000000e-02f * calc_E0_d));                                                               \
    }
#define IFNUMBER_7(name)                                                                                               \
    if ((V_old_ <= 0.000000000000000e+00f)) {                                                                          \
         (name) = (((-8.750000000000000e-04f) * V_old_) + 1.200000000000000e-02f);                                     \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = 1.200000000000000e-02f;                                                                               \
    }
#define IFNUMBER_8(name)                                                                                               \
    if ((calc_dVdt > 1.500000000000000e+05f)) {                                                                        \
         (name) = ((1.000000000000000e+05f * (1.000000000000000e+00f - APtrack_old_)) -                                \
                   (5.000000000000000e+02f * APtrack_old_));                                                           \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((-5.000000000000000e+02f) * APtrack_old_);                                                           \
    }
#define IFNUMBER_9(name)                                                                                               \
    if (((APtrack_old_ < 2.000000000000000e-01f) && (APtrack_old_ > 1.800000000000000e-01f))) {                        \
         (name) = ((1.000000000000000e+05f * (1.000000000000000e+00f - APtrack2_old_)) -                               \
                   (5.000000000000000e+02f * APtrack2_old_));                                                          \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((-5.000000000000000e+02f) * APtrack2_old_);                                                          \
    }
#define IFNUMBER_10(name)                                                                                              \
    if (((APtrack_old_ < 2.000000000000000e-01f) && (APtrack_old_ > 1.800000000000000e-01f))) {                        \
         (name) = ((1.000000000000000e+05f * (1.000000000000000e+00f - APtrack3_old_)) -                               \
                   (5.000000000000000e+02f * APtrack3_old_));                                                          \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((-1.000000000000000e+01f) * APtrack3_old_);                                                          \
    }
#define IFNUMBER_11(name)                                                                                              \
    if ((APtrack_old_ > 2.000000000000000e-01f)) {                                                                     \
         (name) = (((-A_cap) * (((calc_i_CaCa + calc_i_Ca_T) - calc_i_NaCa) + calc_i_p_Ca + calc_i_Ca_b)) /            \
                   (2.000000000000000e+00f * calc_V_myo * F));                                                         \
    }                                                                                                                  \
    else if (((APtrack2_old_ > 1.000000000000000e-02f) && (APtrack_old_ <= 2.000000000000000e-01f))) {                 \
         (name) = 0.000000000000000e+00f;                                                                              \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((-5.000000000000000e+02f) * Cainfluxtrack_old_);                                                     \
    }
#define IFNUMBER_12(name)                                                                                              \
    if ((((1.000000000000000e+00f / (1.000000000000000e+00f + (K_mCSQN / Ca_JSR_old_))) > CSQNthresh) &&               \
         (OVRLDtrack3_old_ < 3.700000000000000e-01f) && (APtrack3_old_ < 3.700000000000000e-01f))) {                   \
         (name) = (5.000000000000000e+04f * (1.000000000000000e+00f - OVRLDtrack_old_));                               \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((-5.000000000000000e+02f) * OVRLDtrack_old_);                                                        \
    }
#define IFNUMBER_13(name)                                                                                              \
    if (((OVRLDtrack_old_ > Logicthresh) && (OVRLDtrack2_old_ < Logicthresh))) {                                       \
         (name) = (5.000000000000000e+04f * (1.000000000000000e+00f - OVRLDtrack2_old_));                              \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((-5.000000000000000e+02f) * OVRLDtrack2_old_);                                                       \
    }
#define IFNUMBER_14(name)                                                                                              \
    if (((OVRLDtrack_old_ > Logicthresh) && (OVRLDtrack3_old_ < Logicthresh))) {                                       \
         (name) = (5.000000000000000e+04f * (1.000000000000000e+00f - OVRLDtrack3_old_));                              \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = ((-1.000000000000000e+01f) * OVRLDtrack3_old_);                                                       \
    }
#define IFNUMBER_15(name)                                                                                              \
    if ((Cainfluxtrack_old_ > delta_Ca_ith)) {                                                                         \
         (name) =                                                                                                      \
            (((G_rel_max * (Cainfluxtrack_old_ - delta_Ca_ith)) / ((K_mrel + Cainfluxtrack_old_) - delta_Ca_ith)) *    \
             (1.000000000000000e+00f - APtrack2_old_) * APtrack2_old_);                                                \
    }                                                                                                                  \
    else if (((Cainfluxtrack_old_ <= delta_Ca_ith) && (OVRLDtrack2_old_ > 0.000000000000000e+00f))) {                  \
         (name) = (G_rel_overload * (1.000000000000000e+00f - OVRLDtrack2_old_) * OVRLDtrack2_old_);                   \
    }                                                                                                                  \
    else {                                                                                                             \
        (name) = 0.000000000000000e+00f;                                                                               \
    }

void RHS_cpu (const real *sv, real *rDY_, real stim_current) {

    // State variables
    const real V_old_ = sv[0];
    const real m_old_ = sv[1];
    const real h_old_ = sv[2];
    const real j_old_ = sv[3];
    const real d_old_ = sv[4];
    const real f_old_ = sv[5];
    const real b_old_ = sv[6];
    const real g_old_ = sv[7];
    const real xr_old_ = sv[8];
    const real xs1_old_ = sv[9];
    const real xs2_old_ = sv[10];
    const real zdv_old_ = sv[11];
    const real ydv_old_ = sv[12];
    const real APtrack_old_ = sv[13];
    const real APtrack2_old_ = sv[14];
    const real APtrack3_old_ = sv[15];
    const real Cainfluxtrack_old_ = sv[16];
    const real OVRLDtrack_old_ = sv[17];
    const real OVRLDtrack2_old_ = sv[18];
    const real OVRLDtrack3_old_ = sv[19];
    const real Ca_JSR_old_ = sv[20];
    const real Ca_NSR_old_ = sv[21];
    const real Cai_old_ = sv[22];
    const real Nai_old_ = sv[23];
    const real Ki_old_ = sv[24];

    // Parameters
    const real Cm = 1.000000000000000e-03f;
    const real R = 8.314000000000000e+03f;
    const real T = 3.100000000000000e+02f;
    const real F = 9.648500000000000e+04f;
    const real Nao = 1.320000000000000e+02f;
    const real g_Na = 1.600000000000000e+01f;
    const real delta_m = 1.000000000000000e-05f;
    const real P_Ca = 5.400000000000000e-04f;
    const real gamma_Cai = 1.000000000000000e+00f;
    const real gamma_Cao = 3.410000000000000e-01f;
    const real Cao = 1.800000000000000e+00f;
    const real P_Na = 6.750000000000000e-07f;
    const real gamma_Nai = 7.500000000000000e-01f;
    const real gamma_Nao = 7.500000000000000e-01f;
    const real P_K = 1.930000000000000e-07f;
    const real gamma_Ki = 7.500000000000000e-01f;
    const real gamma_Ko = 7.500000000000000e-01f;
    const real Ko = 4.500000000000000e+00f;
    const real Km_Ca = 5.999999999999999e-04f;
    const real g_CaT = 5.000000000000000e-02f;
    const real PNaK = 1.833000000000000e-02f;
    const real g_Kp = 5.520000000000000e-03f;
    const real kdKNa = 6.600000000000000e+01f;
    const real nKNa = 2.800000000000000e+00f;
    const real i_K_ATP_on = 1.000000000000000e+00f;
    const real nicholsarea = 5.000000000000000e-05f;
    const real ATPi = 3.000000000000000e+00f;
    const real kATP = 2.500000000000000e-04f;
    const real hATP = 2.000000000000000e+00f;
    const real nATP = 2.400000000000000e-01f;
    const real I_pCa = 1.150000000000000e+00f;
    const real K_mpCa = 5.000000000000000e-04f;
    const real g_Nab = 4.000000000000000e-03f;
    const real g_Cab = 3.016000000000000e-03f;
    const real I_NaK = 2.250000000000000e+00f;
    const real K_mNai = 1.000000000000000e+01f;
    const real K_mKo = 1.500000000000000e+00f;
    const real K_m_ns_Ca = 1.200000000000000e-03f;
    const real c1 = 2.500000000000000e-04f;
    const real gamma = 1.500000000000000e-01f;
    const real c2 = 1.000000000000000e-04f;
    const real A_cap = 1.434000000000000e-07f;
    const real K_mCSQN = 8.000000000000000e-01f;
    const real CSQNthresh = 7.000000000000000e-01f;
    const real Logicthresh = 9.800000000000000e-01f;
    const real G_rel_max = 6.000000000000000e+04f;
    const real delta_Ca_ith = 1.800000000000000e-04f;
    const real K_mrel = 8.000000000000000e-04f;
    const real G_rel_overload = 4.000000000000000e+03f;
    const real I_up = 8.750000000000000e+00f;
    const real K_mup = 9.200000000000000e-04f;
    const real Ca_NSR_max = 1.500000000000000e+01f;
    const real tau_tr = 1.800000000000000e-01f;
    const real CSQN_max = 1.000000000000000e+01f;
    const real CMDN_max = 5.000000000000000e-02f;
    const real K_mCMDN = 2.380000000000000e-03f;
    const real Tn_max = 7.000000000000001e-02f;
    const real K_mTn = 5.000000000000000e-04f;
    const real preplength = 1.000000000000000e-03f;
    const real radius = 1.100000000000000e-04f;

    real calc_E_Na = (((R * T) / F) * logf ((Nao / Nai_old_)));                                // 3
    real calc_E0_m = (V_old_ + 4.713000000000000e+01f);                                        // 5
    real calc_beta_m = (8.000000000000000e+01f * expf (((-V_old_) / 1.100000000000000e+01f))); // 7
    real calc_alpha_h = 0.0f;
    IFNUMBER_2 (calc_alpha_h); // 9
    real calc_beta_h = 0.0f;
    IFNUMBER_3 (calc_beta_h); // 10
    real calc_alpha_j = 0.0f;
    IFNUMBER_4 (calc_alpha_j); // 12
    real calc_beta_j = 0.0f;
    IFNUMBER_5 (calc_beta_j); // 13
    real calc_I_CaCa =
        ((((P_Ca * powf (2.000000000000000e+00f, 2.000000000000000e+00f) * V_old_ * powf (F, 2.000000000000000e+00f)) /
           (R * T)) *
          ((gamma_Cai * Cai_old_ * expf (((2.000000000000000e+00f * V_old_ * F) / (R * T)))) - (gamma_Cao * Cao))) /
         (expf (((2.000000000000000e+00f * V_old_ * F) / (R * T))) - 1.000000000000000e+00f)); // 15
    real calc_I_CaNa =
        ((((P_Na * powf (1.000000000000000e+00f, 2.000000000000000e+00f) * V_old_ * powf (F, 2.000000000000000e+00f)) /
           (R * T)) *
          ((gamma_Nai * Nai_old_ * expf (((1.000000000000000e+00f * V_old_ * F) / (R * T)))) - (gamma_Nao * Nao))) /
         (expf (((1.000000000000000e+00f * V_old_ * F) / (R * T))) - 1.000000000000000e+00f)); // 16
    real calc_I_CaK =
        ((((P_K * powf (1.000000000000000e+00f, 2.000000000000000e+00f) * V_old_ * powf (F, 2.000000000000000e+00f)) /
           (R * T)) *
          ((gamma_Ki * Ki_old_ * expf (((1.000000000000000e+00f * V_old_ * F) / (R * T)))) - (gamma_Ko * Ko))) /
         (expf (((1.000000000000000e+00f * V_old_ * F) / (R * T))) - 1.000000000000000e+00f)); // 17
    real calc_E0_d = (V_old_ + 1.000000000000000e+01f);                                        // 22
    real calc_f_infinity =
        ((1.000000000000000e+00f /
          (1.000000000000000e+00f + expf (((V_old_ + 3.200000000000000e+01f) / 8.000000000000000e+00f)))) +
         (6.000000000000000e-01f /
          (1.000000000000000e+00f + expf (((5.000000000000000e+01f - V_old_) / 2.000000000000000e+01f))))); // 28
    real calc_tau_f =
        (1.000000000000000e-03f /
         ((1.970000000000000e-02f *
           expf ((-powf ((3.370000000000000e-02f * (V_old_ + 1.000000000000000e+01f)), 2.000000000000000e+00f)))) +
          2.000000000000000e-02f));                                                            // 29
    real calc_f_Ca = (1.000000000000000e+00f / (1.000000000000000e+00f + (Cai_old_ / Km_Ca))); // 33
    real calc_b_inf =
        (1.000000000000000e+00f /
         (1.000000000000000e+00f + expf (((-(V_old_ + 1.400000000000000e+01f)) / 1.080000000000000e+01f)))); // 35
    real calc_tau_b =
        (3.700000000000000e-03f +
         (6.100000000000000e-03f /
          (1.000000000000000e+00f + expf (((V_old_ + 2.500000000000000e+01f) / 4.500000000000000e+00f))))); // 36
    real calc_g_inf =
        (1.000000000000000e+00f /
         (1.000000000000000e+00f + expf (((V_old_ + 6.000000000000000e+01f) / 5.600000000000000e+00f)))); // 38
    real calc_tau_g = 0.0f;
    IFNUMBER_7 (calc_tau_g);                                                                     // 39
    real calc_g_Kr = (2.614000000000000e-02f * powf ((Ko / 5.400000000000000e+00f), 1.0 / 2.0)); // 41
    real calc_Rect =
        (1.000000000000000e+00f /
         (1.000000000000000e+00f + expf (((V_old_ + 9.000000000000000e+00f) / 2.240000000000000e+01f)))); // 42
    real calc_xr_infinity =
        (1.000000000000000e+00f /
         (1.000000000000000e+00f + expf (((-(V_old_ + 2.150000000000000e+01f)) / 7.500000000000000e+00f)))); // 44
    real calc_tau_xr =
        (1.000000000000000e-03f /
         (((1.380000000000000e-03f * (V_old_ + 1.420000000000000e+01f)) /
           (1.000000000000000e+00f - expf (((-1.230000000000000e-01f) * (V_old_ + 1.420000000000000e+01f))))) +
          ((6.100000000000000e-04f * (V_old_ + 3.890000000000000e+01f)) /
           (expf ((1.450000000000000e-01f * (V_old_ + 3.890000000000000e+01f))) - 1.000000000000000e+00f)))); // 45
    real calc_E_Ks = (((R * T) / F) * logf (((Ko + (PNaK * Nao)) / (Ki_old_ + (PNaK * Nai_old_)))));          // 47
    real calc_g_Ks = (4.330000000000000e-01f *
                      (1.000000000000000e+00f +
                       (6.000000000000000e-01f / (1.000000000000000e+00f + powf ((3.800000000000000e-05f / Cai_old_),
                                                                                 1.400000000000000e+00f))))); // 48
    real calc_xs1_infinity =
        (1.000000000000000e+00f /
         (1.000000000000000e+00f + expf (((-(V_old_ - 1.500000000000000e+00f)) / 1.670000000000000e+01f)))); // 50
    real calc_tau_xs1 =
        (1.000000000000000e-03f /
         (((7.190000000000000e-05f * (V_old_ + 3.000000000000000e+01f)) /
           (1.000000000000000e+00f - expf (((-1.480000000000000e-01f) * (V_old_ + 3.000000000000000e+01f))))) +
          ((1.310000000000000e-04f * (V_old_ + 3.000000000000000e+01f)) /
           (expf ((6.870000000000000e-02f * (V_old_ + 3.000000000000000e+01f))) - 1.000000000000000e+00f)))); // 51
    real calc_xs2_infinity =
        (1.000000000000000e+00f /
         (1.000000000000000e+00f + expf (((-(V_old_ - 1.500000000000000e+00f)) / 1.670000000000000e+01f)))); // 53
    real calc_tau_xs2 =
        ((4.000000000000000e+00f * 1.000000000000000e-03f) /
         (((7.190000000000000e-05f * (V_old_ + 3.000000000000000e+01f)) /
           (1.000000000000000e+00f - expf (((-1.480000000000000e-01f) * (V_old_ + 3.000000000000000e+01f))))) +
          ((1.310000000000000e-04f * (V_old_ + 3.000000000000000e+01f)) /
           (expf ((6.870000000000000e-02f * (V_old_ + 3.000000000000000e+01f))) - 1.000000000000000e+00f)))); // 54
    real calc_g_K1 = (7.500000000000000e-01f * powf ((Ko / 5.400000000000000e+00f), 1.0 / 2.0));              // 56
    real calc_E_K = (((R * T) / F) * logf ((Ko / Ki_old_)));                                                  // 57
    real calc_Kp =
        (1.000000000000000e+00f /
         (1.000000000000000e+00f + expf (((7.488000000000000e+00f - V_old_) / 5.980000000000000e+00f))));   // 62
    real calc_g_K_Na = (0.000000000000000e+00f * 1.284800000000000e-01f);                                   // 64
    real calc_pona = (8.500000000000000e-01f / (1.000000000000000e+00f + powf ((kdKNa / Nai_old_), nKNa))); // 65
    real calc_pov =
        (8.000000000000000e-01f -
         (6.500000000000000e-01f /
          (1.000000000000000e+00f + expf (((V_old_ + 1.250000000000000e+02f) / 1.500000000000000e+01f))))); // 66
    real calc_g_K_ATP = ((i_K_ATP_on * 1.930000000000000e-04f) / nicholsarea);                              // 68
    real calc_pATP = (1.000000000000000e+00f / (1.000000000000000e+00f + powf ((ATPi / kATP), hATP)));      // 69
    real calc_g_to = (0.000000000000000e+00f * 5.000000000000000e-01f);                                     // 72
    real calc_rvdv = expf ((V_old_ / 1.000000000000000e+02f));                                              // 73
    real calc_alpha_zdv =
        ((1.000000000000000e+04f * expf (((V_old_ - 4.000000000000000e+01f) / 2.500000000000000e+01f))) /
         (1.000000000000000e+00f + expf (((V_old_ - 4.000000000000000e+01f) / 2.500000000000000e+01f)))); // 75
    real calc_beta_zdv =
        ((1.000000000000000e+04f * expf (((-(V_old_ + 9.000000000000000e+01f)) / 2.500000000000000e+01f))) /
         (1.000000000000000e+00f + expf (((-(V_old_ + 9.000000000000000e+01f)) / 2.500000000000000e+01f)))); // 76
    real calc_alpha_ydv =
        (1.500000000000000e+01f /
         (1.000000000000000e+00f + expf (((V_old_ + 6.000000000000000e+01f) / 5.000000000000000e+00f)))); // 80
    real calc_beta_ydv =
        ((1.000000000000000e+02f * expf (((V_old_ + 2.500000000000000e+01f) / 5.000000000000000e+00f))) /
         (1.000000000000000e+00f + expf (((V_old_ + 2.500000000000000e+01f) / 5.000000000000000e+00f)))); // 81
    real calc_i_p_Ca = ((I_pCa * Cai_old_) / (K_mpCa + Cai_old_));                                        // 85
    real calc_E_Ca = (((R * T) / (2.000000000000000e+00f * F)) * logf ((Cao / Cai_old_)));                // 87
    real calc_sigma = ((1.000000000000000e+00f / 7.000000000000000e+00f) *
                       (expf ((Nao / 6.730000000000000e+01f)) - 1.000000000000000e+00f)); // 89
    real calc_P_ns_Ca = (0.000000000000000e+00f * 1.750000000000000e-07f);                // 92
    real calc_i_NaCa =
        ((c1 * expf ((((gamma - 1.000000000000000e+00f) * V_old_ * F) / (R * T))) *
          ((expf (((V_old_ * F) / (R * T))) * powf (Nai_old_, 3.000000000000000e+00f) * Cao) -
           (powf (Nao, 3.000000000000000e+00f) * Cai_old_))) /
         (1.000000000000000e+00f + (c2 * expf ((((gamma - 1.000000000000000e+00f) * V_old_ * F) / (R * T))) *
                                    ((expf (((V_old_ * F) / (R * T))) * powf (Nai_old_, 3.000000000000000e+00f) * Cao) +
                                     (powf (Nao, 3.000000000000000e+00f) * Cai_old_))))); // 98
    real calc_G_rel = 0.0f;
    IFNUMBER_15 (calc_G_rel);                                                                                 // 108
    real calc_i_up = ((I_up * Cai_old_) / (Cai_old_ + K_mup));                                                // 110
    real calc_K_leak = (I_up / Ca_NSR_max);                                                                   // 111
    real calc_i_tr = ((Ca_NSR_old_ - Ca_JSR_old_) / tau_tr);                                                  // 113
    real calc_volume = (3.141592653589793116f * preplength * powf (radius, 2.000000000000000e+00f));          // 117
    real calc_i_Na = (g_Na * powf (m_old_, 3.000000000000000e+00f) * h_old_ * j_old_ * (V_old_ - calc_E_Na)); // 4
    real calc_alpha_m = 0.0f;
    IFNUMBER_1 (calc_alpha_m); // 6
    real calc_d_infinity =
        (1.000000000000000e+00f / (1.000000000000000e+00f + expf (((-calc_E0_d) / 6.240000000000000e+00f)))); // 23
    real calc_alpha_f = (calc_f_infinity / calc_tau_f);                                                       // 30
    real calc_beta_f = ((1.000000000000000e+00f - calc_f_infinity) / calc_tau_f);                             // 31
    real calc_i_Ks = (calc_g_Ks * xs1_old_ * xs2_old_ * (V_old_ - calc_E_Ks));                                // 49
    real calc_alpha_K1 = (1.020000000000000e+03f /
                          (1.000000000000000e+00f +
                           expf ((2.385000000000000e-01f * ((V_old_ - calc_E_K) - 5.921500000000000e+01f))))); // 59
    real calc_beta_K1 =
        ((1.000000000000000e+03f *
          ((4.912400000000000e-01f * expf ((8.032000000000000e-02f * ((V_old_ - calc_E_K) + 5.476000000000000e+00f)))) +
           expf ((6.175000000000000e-02f * ((V_old_ - calc_E_K) - 5.943099999999999e+02f))))) /
         (1.000000000000000e+00f +
          expf (((-5.143000000000000e-01f) * ((V_old_ - calc_E_K) + 4.753000000000000e+00f)))));   // 60
    real calc_i_Kp = (g_Kp * calc_Kp * (V_old_ - calc_E_K));                                       // 63
    real calc_i_K_Na = (calc_g_K_Na * calc_pona * calc_pov * (V_old_ - calc_E_K));                 // 67
    real calc_GKbaraATP = (calc_g_K_ATP * calc_pATP * powf ((Ko / 4.000000000000000e+00f), nATP)); // 70
    real calc_i_to =
        (calc_g_to * powf (zdv_old_, 3.000000000000000e+00f) * ydv_old_ * calc_rvdv * (V_old_ - calc_E_K)); // 74
    real calc_tau_zdv = (1.000000000000000e+00f / (calc_alpha_zdv + calc_beta_zdv));                        // 77
    real calc_zdv_ss = (calc_alpha_zdv / (calc_alpha_zdv + calc_beta_zdv));                                 // 78
    real calc_tau_ydv = (1.000000000000000e+00f / (calc_alpha_ydv + calc_beta_ydv));                        // 82
    real calc_ydv_ss = (calc_alpha_ydv / (calc_alpha_ydv + calc_beta_ydv));                                 // 83
    real calc_i_Na_b = (g_Nab * (V_old_ - calc_E_Na));                                                      // 86
    real calc_i_Ca_b = (g_Cab * (V_old_ - calc_E_Ca));                                                      // 88
    real calc_f_NaK = (1.000000000000000e+00f /
                       (1.000000000000000e+00f +
                        (1.245000000000000e-01f * expf ((((-1.000000000000000e-01f) * V_old_ * F) / (R * T)))) +
                        (3.650000000000000e-02f * calc_sigma * expf ((((-V_old_) * F) / (R * T)))))); // 90
    real calc_I_ns_Na =
        ((((calc_P_ns_Ca * powf (1.000000000000000e+00f, 2.000000000000000e+00f) * V_old_ *
            powf (F, 2.000000000000000e+00f)) /
           (R * T)) *
          ((gamma_Nai * Nai_old_ * expf (((1.000000000000000e+00f * V_old_ * F) / (R * T)))) - (gamma_Nao * Nao))) /
         (expf (((1.000000000000000e+00f * V_old_ * F) / (R * T))) - 1.000000000000000e+00f)); // 93
    real calc_I_ns_K =
        ((((calc_P_ns_Ca * powf (1.000000000000000e+00f, 2.000000000000000e+00f) * V_old_ *
            powf (F, 2.000000000000000e+00f)) /
           (R * T)) *
          ((gamma_Ki * Ki_old_ * expf (((1.000000000000000e+00f * V_old_ * F) / (R * T)))) - (gamma_Ko * Ko))) /
         (expf (((1.000000000000000e+00f * V_old_ * F) / (R * T))) - 1.000000000000000e+00f)); // 94
    real calc_i_rel = (calc_G_rel * (Ca_JSR_old_ - Cai_old_));                                 // 109
    real calc_i_leak = (calc_K_leak * Ca_NSR_old_);                                            // 112
    real calc_V_myo = (6.800000000000000e-01f * calc_volume);                                  // 118
    real calc_i_CaCa = (d_old_ * f_old_ * calc_f_Ca * calc_I_CaCa);                            // 18
    real calc_i_CaNa = (d_old_ * f_old_ * calc_f_Ca * calc_I_CaNa);                            // 19
    real calc_i_CaK = (d_old_ * f_old_ * calc_f_Ca * calc_I_CaK);                              // 20
    real calc_tau_d = 0.0f;
    IFNUMBER_6 (calc_tau_d);                                                      // 24
    real calc_i_Ca_T = (g_CaT * b_old_ * b_old_ * g_old_ * (V_old_ - calc_E_Ca)); // 34
    real calc_i_Kr = (calc_g_Kr * xr_old_ * calc_Rect * (V_old_ - calc_E_K));     // 43
    real calc_K1_infinity = (calc_alpha_K1 / (calc_alpha_K1 + calc_beta_K1));     // 61
    real calc_i_K_ATP = (calc_GKbaraATP * (V_old_ - calc_E_K));                   // 71
    real calc_i_NaK = ((((I_NaK * calc_f_NaK * 1.000000000000000e+00f) /
                         (1.000000000000000e+00f + powf ((K_mNai / Nai_old_), 2.000000000000000e+00f))) *
                        Ko) /
                       (Ko + K_mKo)); // 91
    real calc_i_ns_Na = ((calc_I_ns_Na * 1.000000000000000e+00f) /
                         (1.000000000000000e+00f + powf ((K_m_ns_Ca / Cai_old_), 3.000000000000000e+00f))); // 95
    real calc_i_ns_K = ((calc_I_ns_K * 1.000000000000000e+00f) /
                        (1.000000000000000e+00f + powf ((K_m_ns_Ca / Cai_old_), 3.000000000000000e+00f))); // 96
    real calc_i_Ca_L = (calc_i_CaCa + calc_i_CaK + calc_i_CaNa);                                           // 21
    real calc_alpha_d = (calc_d_infinity / calc_tau_d);                                                    // 25
    real calc_beta_d = ((1.000000000000000e+00f - calc_d_infinity) / calc_tau_d);                          // 26
    real calc_i_ns_Ca = (calc_i_ns_Na + calc_i_ns_K);                                                      // 97
    real calc_V_JSR = ((4.800000000000000e-03f / 6.800000000000000e-01f) * calc_V_myo);                    // 99
    real calc_V_NSR = ((5.520000000000000e-02f / 6.800000000000000e-01f) * calc_V_myo);                    // 100
    real calc_i_K1 = (calc_g_K1 * calc_K1_infinity * (V_old_ - calc_E_K));                                 // 58
    real calc_dVdt = (((-1.000000000000000e+00f) / Cm) *
                      (calc_i_Na + calc_i_Ca_L + calc_i_Ca_T + calc_i_Kr + calc_i_Ks + calc_i_K_Na + calc_i_K_ATP +
                       calc_i_to + calc_i_K1 + calc_i_Kp + calc_i_NaCa + calc_i_p_Ca + calc_i_Na_b + calc_i_Ca_b +
                       calc_i_NaK + calc_i_ns_Ca + stim_current));

    rDY_[0] = calc_dVdt;
    rDY_[1] = ((calc_alpha_m * (1.000000000000000e+00f - m_old_)) - (calc_beta_m * m_old_));
    rDY_[2] = ((calc_alpha_h * (1.000000000000000e+00f - h_old_)) - (calc_beta_h * h_old_));
    rDY_[3] = ((calc_alpha_j * (1.000000000000000e+00f - j_old_)) - (calc_beta_j * j_old_));
    rDY_[4] = ((calc_alpha_d * (1.000000000000000e+00f - d_old_)) - (calc_beta_d * d_old_));
    rDY_[5] = ((calc_alpha_f * (1.000000000000000e+00f - f_old_)) - (calc_beta_f * f_old_));
    rDY_[6] = ((calc_b_inf - b_old_) / calc_tau_b);
    rDY_[7] = ((calc_g_inf - g_old_) / calc_tau_g);
    rDY_[8] = ((calc_xr_infinity - xr_old_) / calc_tau_xr);
    rDY_[9] = ((calc_xs1_infinity - xs1_old_) / calc_tau_xs1);
    rDY_[10] = ((calc_xs2_infinity - xs2_old_) / calc_tau_xs2);
    rDY_[11] = ((calc_zdv_ss - zdv_old_) / calc_tau_zdv);
    rDY_[12] = ((calc_ydv_ss - ydv_old_) / calc_tau_ydv);
    rDY_[13] = 0.0f;
    IFNUMBER_8 (rDY_[13]);
    rDY_[14] = 0.0f;
    IFNUMBER_9 (rDY_[14]);
    rDY_[15] = 0.0f;
    IFNUMBER_10 (rDY_[15]);
    rDY_[16] = 0.0f;
    IFNUMBER_11 (rDY_[16]);
    rDY_[17] = 0.0f;
    IFNUMBER_12 (rDY_[17]);
    rDY_[18] = 0.0f;
    IFNUMBER_13 (rDY_[18]);
    rDY_[19] = 0.0f;
    IFNUMBER_14 (rDY_[19]);
    rDY_[20] =
        ((1.000000000000000e+00f /
          (1.000000000000000e+00f + ((CSQN_max * K_mCSQN) / powf ((K_mCSQN + Ca_JSR_old_), 2.000000000000000e+00f)))) *
         (calc_i_tr - calc_i_rel));
    rDY_[21] = (((((-calc_i_tr) * calc_V_JSR) / calc_V_NSR) - calc_i_leak) + calc_i_up);
    rDY_[22] =
        ((1.000000000000000e+00f /
          (1.000000000000000e+00f + ((CMDN_max * K_mCMDN) / powf ((K_mCMDN + Cai_old_), 2.000000000000000e+00f)) +
           ((Tn_max * K_mTn) / powf ((K_mTn + Cai_old_), 2.000000000000000e+00f)))) *
         ((((-A_cap) * (((calc_i_CaCa + calc_i_Ca_T) - calc_i_NaCa) + calc_i_p_Ca + calc_i_Ca_b)) /
           (2.000000000000000e+00f * calc_V_myo * F)) +
          ((calc_i_rel * calc_V_JSR) / calc_V_myo) + (((calc_i_leak - calc_i_up) * calc_V_NSR) / calc_V_myo)));
    rDY_[23] = (((-(calc_i_Na + calc_i_CaNa + calc_i_Na_b + calc_i_ns_Na + (calc_i_NaCa * 3.000000000000000e+00f) +
                    (calc_i_NaK * 3.000000000000000e+00f))) *
                 A_cap) /
                (calc_V_myo * F));
    rDY_[24] = (((-(calc_i_CaK + calc_i_Kr + calc_i_Ks + calc_i_K1 + calc_i_Kp + calc_i_K_Na + calc_i_K_ATP +
                    calc_i_to + calc_i_ns_K + ((-calc_i_NaK) * 2.000000000000000e+00f))) *
                 A_cap) /
                (calc_V_myo * F));
}
