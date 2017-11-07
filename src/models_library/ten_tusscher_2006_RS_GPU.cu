#include <stddef.h>
#include <stdint.h>
#include "model_gpu_utils.h"

#include "ten_tusscher_2006.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    print_to_stdout_and_file("Using ten Tusscher 2006 GPU model\n");

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}


extern "C" SOLVE_MODEL_ODES_GPU(solve_model_odes_gpu) {


    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;


    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));


    //the array cells to solve is passed when we are using and adapative mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }
    solve_gpu <<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes)
{
    // Thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if(threadID < num_volumes) {

        *((real*)((char*)sv + pitch * 0) + threadID) = -85.23f;   // V;       millivolt
        *((real*)((char*)sv + pitch * 1) + threadID) = 0.00621;  // Xr1;     dimensionless
        *((real*)((char*)sv + pitch * 2) + threadID) = 0.4712;   // Xr2;     dimensionless
        *((real*)((char*)sv + pitch * 3) + threadID) = 0.0095;   // Xs;      dimensionless
        *((real*)((char*)sv + pitch * 4) + threadID) = 0.00172;  // m;       dimensionless
        *((real*)((char*)sv + pitch * 5) + threadID) = 0.7444;   // h;       dimensionless
        *((real*)((char*)sv + pitch * 6) + threadID) = 0.7045;   // j;       dimensionless
        *((real*)((char*)sv + pitch * 7) + threadID) = 3.373e-5; // d;       dimensionless
        *((real*)((char*)sv + pitch * 8) + threadID) = 0.7888;   // f;       dimensionless
        *((real*)((char*)sv + pitch * 9) + threadID) = 0.9755;   // f2;      dimensionless
        *((real*)((char*)sv + pitch * 10) + threadID) = 0.9953;   // fCass;   dimensionless
        *((real*)((char*)sv + pitch * 11) + threadID) = 0.999998; // s;       dimensionless
        *((real*)((char*)sv + pitch * 12) + threadID) = 2.42e-8;  // r;       dimensionless
        *((real*)((char*)sv + pitch * 13) + threadID) = 0.000126; // Ca_i;    millimolar
        *((real*)((char*)sv + pitch * 14) + threadID) = 3.64;     // Ca_SR;   millimolar
        *((real*)((char*)sv + pitch * 15) + threadID) = 0.00036;  // Ca_ss;   millimolar
        *((real*)((char*)sv + pitch * 16) + threadID) = 0.9073;   // R_prime; dimensionless
        *((real*)((char*)sv + pitch * 17) + threadID) = 8.604;    // Na_i;    millimolar
        *((real*)((char*)sv + pitch * 18) + threadID) = 136.89;   // K_i;     millimolar

    }
}


// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt);

            *((real*)((char*)sv) + sv_id) = dt*rDY[0] + *((real*)((char*)sv) + sv_id);

            for(int i = 1; i < 13; i++) {
                *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }

            for(int i = 13; i < 19; i++) {
                *((real *) ((char *) sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *) ((char *) sv + pitch * i) + sv_id);
            }
            
        }

    }
}


inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt) {

    // State variables
    const real V = *((real*)((char*)sv_ + pitch * 0) + threadID_);      // Membrane variable
    const real Xr1 = *((real*)((char*)sv_ + pitch * 1) + threadID_);    // Rapid time dependent potassium current Xr1
    const real Xr2 = *((real*)((char*)sv_ + pitch * 2) + threadID_);    // Rapid time dependent potassium current Xr2
    const real Xs = *((real*)((char*)sv_ + pitch * 3) + threadID_);     // Slow time dependent potassium current Xs
    const real m = *((real*)((char*)sv_ + pitch * 4) + threadID_);      // Fast sodium current m
    const real h = *((real*)((char*)sv_ + pitch * 5) + threadID_);      // Fast sodium current h gate
    const real j = *((real*)((char*)sv_ + pitch * 6) + threadID_);      // Fast sodium current j gate
    const real d = *((real*)((char*)sv_ + pitch * 7) + threadID_);      // L type Ca current d gate
    const real f = *((real*)((char*)sv_ + pitch * 8) + threadID_);;      // var_L_type_Ca_current_f_gate__f
    const real f2 = *((real*)((char*)sv_ + pitch * 9) + threadID_);     // var_L_type_Ca_current_f2_gate__f2
    const real fCass = *((real*)((char*)sv_ + pitch * 10) + threadID_); // L_type_Ca_current__fCass
    const real s = *((real*)((char*)sv_ + pitch * 11) + threadID_);     // gating s
    const real r = *((real*)((char*)sv_ + pitch * 12) + threadID_);     // gating r
    const real Ca_i = *((real*)((char*)sv_ + pitch * 13) + threadID_);  // calcium_dynamics__Ca_i
    const real Ca_SR = *((real*)((char*)sv_ + pitch * 14) + threadID_);
    const real Ca_ss = *((real*)((char*)sv_ + pitch * 15) + threadID_);
    const real R_prime = *((real*)((char*)sv_ + pitch * 16) + threadID_);
    const real Na_i = *((real*)((char*)sv_ + pitch * 17) + threadID_); // var_sodium_dynamics__Na_i
    const real K_i = *((real*)((char*)sv_ + pitch * 18) + threadID_);  // var_potassium_dynamics__K_i

    // Some constants
    const real R   = 8314.472;
    const real T   = 310.0;
    const real F   = 96485.3415f;
    const real Cm  = 0.185;
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
    //real var_membrane__i_Stim = ((time>=stim_start)&&(time<=stim_start+stim_dur)) ? stim_current: 0.0f;
    real var_membrane__i_Stim =  stim_current;

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

/*extern "C" void update_gpu_after_refinement(real *sv, uint32_t *cells, size_t number_of_cells, int neq) {

*//*    // execution configuration
    const int GRID  = ((int)number_of_cells + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = number_of_cells*sizeof(uint32_t);

    uint32_t *cells_d;
    check_cuda_error(cudaMalloc((void **) &cells_d, size));
    check_cuda_error(cudaMemcpy(cells_d, cells, size, cudaMemcpyHostToDevice));

    update_refinement<<<GRID, BLOCK_SIZE>>>(sv, cells_d, number_of_cells, neq);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(cells_d));*//*

    real *sv_src;
    real *sv_dst;

    for (size_t i = 0; i < number_of_cells/8; i++) {

        size_t index_id = i * 8;

        uint32_t index = cells[index_id];
        sv_src = &sv[index];

        for (int j = 1; j < 8; j++) {
            index = cells[index_id + j];
            sv_dst = &sv[index];
            cudaMemcpy2D(sv_dst, pitch_h, sv_src, pitch_h, sizeof(real), (size_t )neq, cudaMemcpyDeviceToDevice);
        }


    }


}*/

/*
__global__ void update_refinement(real *sv, uint32_t *cells, size_t number_of_cells, int neq) {

    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int index_dst;
    int i = 0;
    int index_id = threadID*8;
    int index_src;

    if(index_id < number_of_cells) {

        index_src = cells[index_id];
//        real *src = (real*)malloc(sizeof(real)*neq);
//
//        for(int k = 0; k < neq; k++) {
//            src[i] = *((real*)((char*)sv + pitch * k) + index_dst);
//        }

        for(i = 1; i < 8; i++) {

            index_dst = cells[index_id+i];

            for(int k = 0; k < neq; k++) {
                *((real*)((char*)sv + pitch * k) + index_dst) = *((real*)((char*)sv + pitch * k) + index_src);
            }

        }
    }

}*/
