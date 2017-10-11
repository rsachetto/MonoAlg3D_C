#include <stddef.h>
#include <unitypes.h>
#include "../solvers/constants.h"
#include <stdlib.h>
#include <stdio.h>

/*! \brief Macros/inlines to assist CLion to parse Cuda files (*.cu, *.cuh) */
#ifdef __JETBRAINS_IDE__
#define __CUDACC__ 1
#define __host__
#define __device__
#define __global__
#define __forceinline__
#define __shared__
inline void __syncthreads() {}
inline void __threadfence_block() {}
template<class T> inline T __clz(const T val) { return val; }
struct __cuda_fake_struct { int x; int y; int z; };
extern __cuda_fake_struct blockDim;
extern __cuda_fake_struct threadIdx;
extern __cuda_fake_struct blockIdx;
#include "../../../../../opt/cuda/include/driver_types.h"
#include "../../../../../opt/cuda/include/cuda_runtime.h"
#endif

//TODO: the device and the kernels have to see this value
__device__ size_t pitch;
size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(Real *sv, int num_volumes);

__global__ void solve_gpu(Real dt, Real *sv, Real* stim_currents,
                          uint32_t *cells_to_solve, int num_cells_to_solve,
                          Real stim_start, Real stim_dur, Real time,
                          int num_steps, int neq, void *extra_data);

__global__ void update_refinement(Real *sv, uint32_t *cells, size_t number_of_cells, int neq);

inline __device__ void RHS_gpu(Real *sv_, Real *rDY_, Real stim_current, Real time, Real stim_start, Real stim_dur, int threadID_, Real dt);

void check_err(cudaError_t err) {
    if(err != cudaSuccess) {
        printf("Error code %d\n", err);
        exit(0);
    }

}

extern "C" void set_model_initial_conditions_gpu(Real **sv, uint32_t num_volumes, int neq) {

    printf("Using ten Tusscher GPU model\n");

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(Real);

    // ODE arrays
    check_err(cudaMallocPitch ((void**) &(*sv), &pitch_h, size, neq));
    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_err(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

}

extern "C" void solve_model_ode_gpu(Real dt, Real *sv, Real *stim_currents, uint32_t *cells_to_solve,
                                    uint32_t num_cells_to_solve, Real stim_start, Real stim_dur,
                                    Real time, int num_steps, int neq, void *extra_data) {


    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(Real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    Real *stims_currents_device;
    uint32_t *cells_to_solve_device;
    check_err(cudaMalloc((void**)&stims_currents_device, stim_currents_size));
    check_err(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    check_err(cudaMalloc((void**)&cells_to_solve_device, cells_to_solve_size));
    check_err(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));

    solve_gpu<<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve,
            stim_start, stim_dur, time, num_steps, neq, extra_data);

    check_err(cudaFree(stims_currents_device));
    check_err(cudaFree(cells_to_solve_device));
}


__global__ void kernel_set_model_inital_conditions(Real *sv, int num_volumes)
{
    // Thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if(threadID < num_volumes) {

        *((Real*)((char*)sv + pitch * 0) + threadID) = -85.23f;   // V;       millivolt
        *((Real*)((char*)sv + pitch * 1) + threadID) = 0.00621;  // Xr1;     dimensionless
        *((Real*)((char*)sv + pitch * 2) + threadID) = 0.4712;   // Xr2;     dimensionless
        *((Real*)((char*)sv + pitch * 3) + threadID) = 0.0095;   // Xs;      dimensionless
        *((Real*)((char*)sv + pitch * 4) + threadID) = 0.00172;  // m;       dimensionless
        *((Real*)((char*)sv + pitch * 5) + threadID) = 0.7444;   // h;       dimensionless
        *((Real*)((char*)sv + pitch * 6) + threadID) = 0.7045;   // j;       dimensionless
        *((Real*)((char*)sv + pitch * 7) + threadID) = 3.373e-5; // d;       dimensionless
        *((Real*)((char*)sv + pitch * 8) + threadID) = 0.7888;   // f;       dimensionless
        *((Real*)((char*)sv + pitch * 9) + threadID) = 0.9755;   // f2;      dimensionless
        *((Real*)((char*)sv + pitch * 10) + threadID) = 0.9953;   // fCass;   dimensionless
        *((Real*)((char*)sv + pitch * 11) + threadID) = 0.999998; // s;       dimensionless
        *((Real*)((char*)sv + pitch * 12) + threadID) = 2.42e-8;  // r;       dimensionless
        *((Real*)((char*)sv + pitch * 13) + threadID) = 0.000126; // Ca_i;    millimolar
        *((Real*)((char*)sv + pitch * 14) + threadID) = 3.64;     // Ca_SR;   millimolar
        *((Real*)((char*)sv + pitch * 15) + threadID) = 0.00036;  // Ca_ss;   millimolar
        *((Real*)((char*)sv + pitch * 16) + threadID) = 0.9073;   // R_prime; dimensionless
        *((Real*)((char*)sv + pitch * 17) + threadID) = 8.604;    // Na_i;    millimolar
        *((Real*)((char*)sv + pitch * 18) + threadID) = 136.89;   // K_i;     millimolar

    }
}


// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(Real dt, Real *sv, Real* stim_currents,
                          uint32_t *cells_to_solve, int num_cells_to_solve,
                          Real stim_start, Real stim_dur, Real time,
                          int num_steps, int neq, void *extra_data)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;
    Real t = time;
    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        sv_id = cells_to_solve[threadID];

        for (int i = 0; i < num_steps; ++i) {
            Real *rDY = (Real *)malloc(neq*sizeof(Real));

            RHS_gpu(sv, rDY, stim_currents[threadID], t, stim_start, stim_dur, sv_id, dt);

            *((Real*)((char*)sv) + sv_id) = dt*rDY[0] + *((Real*)((char*)sv) + sv_id);

            for(int j = 1; j < 13; j++) {
                *((Real*)((char*)sv + pitch * j) + sv_id) = rDY[i];
            }

            for(int k = 13; k < 19; k++) {
                *((Real *) ((char *) sv + pitch * k) + sv_id) = dt * rDY[i] + *((Real *) ((char *) sv + pitch * i) + sv_id);
            }
            
            t += dt;
        }

    }
}


inline __device__ void RHS_gpu(Real *sv_, Real *rDY_, Real stim_current, Real time, Real stim_start, Real stim_dur, int threadID_, Real dt) {

    // State variables
    const Real V = *((Real*)((char*)sv_ + pitch * 0) + threadID_);      // Membrane variable
    const Real Xr1 = *((Real*)((char*)sv_ + pitch * 1) + threadID_);    // Rapid time dependent potassium current Xr1
    const Real Xr2 = *((Real*)((char*)sv_ + pitch * 2) + threadID_);    // Rapid time dependent potassium current Xr2
    const Real Xs = *((Real*)((char*)sv_ + pitch * 3) + threadID_);     // Slow time dependent potassium current Xs
    const Real m = *((Real*)((char*)sv_ + pitch * 4) + threadID_);      // Fast sodium current m
    const Real h = *((Real*)((char*)sv_ + pitch * 5) + threadID_);      // Fast sodium current h gate
    const Real j = *((Real*)((char*)sv_ + pitch * 6) + threadID_);      // Fast sodium current j gate
    const Real d = *((Real*)((char*)sv_ + pitch * 7) + threadID_);      // L type Ca current d gate
    const Real f = *((Real*)((char*)sv_ + pitch * 8) + threadID_);;      // var_L_type_Ca_current_f_gate__f
    const Real f2 = *((Real*)((char*)sv_ + pitch * 9) + threadID_);     // var_L_type_Ca_current_f2_gate__f2
    const Real fCass = *((Real*)((char*)sv_ + pitch * 10) + threadID_); // L_type_Ca_current__fCass
    const Real s = *((Real*)((char*)sv_ + pitch * 11) + threadID_);     // gating s
    const Real r = *((Real*)((char*)sv_ + pitch * 12) + threadID_);     // gating r
    const Real Ca_i = *((Real*)((char*)sv_ + pitch * 13) + threadID_);  // calcium_dynamics__Ca_i
    const Real Ca_SR = *((Real*)((char*)sv_ + pitch * 14) + threadID_);
    const Real Ca_ss = *((Real*)((char*)sv_ + pitch * 15) + threadID_);
    const Real R_prime = *((Real*)((char*)sv_ + pitch * 16) + threadID_);
    const Real Na_i = *((Real*)((char*)sv_ + pitch * 17) + threadID_); // var_sodium_dynamics__Na_i
    const Real K_i = *((Real*)((char*)sv_ + pitch * 18) + threadID_);  // var_potassium_dynamics__K_i

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
    Real var_membrane__i_Stim = ((time>=stim_start)&&(time<=stim_start+stim_dur)) ? stim_current: 0.0f;

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

extern "C" void update_gpu_after_refinement(Real *sv, uint32_t *cells, size_t number_of_cells, int neq) {

    // execution configuration
    const int GRID  = ((int)number_of_cells + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = number_of_cells*sizeof(uint32_t);

    uint32_t *cells_d;
    check_err(cudaMalloc((void**)&cells_d, size));
    check_err(cudaMemcpy(cells_d, cells, size, cudaMemcpyHostToDevice));

        update_refinement<<<GRID, BLOCK_SIZE>>>(sv, cells_d, number_of_cells, neq);

    check_err(cudaFree(cells_d));
}

__global__ void update_refinement(Real *sv, uint32_t *cells, size_t number_of_cells, int neq) {

    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int index;
    int i = 0;
    int index_id = threadID*8;

    if(index_id < number_of_cells) {

        index = cells[index_id];
        Real *src = (Real*)malloc(sizeof(Real)*neq);

        for(int k = 0; k < neq; k++) {
            src[i] = *((Real*)((char*)sv + pitch * k) + index);
        }

        for(i = 1; i < 8; i++) {

            index = cells[index_id+i];

            for(int k = 0; k < neq; k++) {
                *((Real*)((char*)sv + pitch * k) + index) = src[i];
            }

        }
    }

}