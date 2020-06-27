#include <stddef.h>
#include <stdint.h>
#include "model_gpu_utils.h"

#include "stewart_aslanidi_noble_2009.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    log_to_stdout_and_file("Using Stewart-Aslanidi-Noble 2009 GPU model\n");

    uint32_t num_volumes = solver->original_num_cells;

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(solver->sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_gpu) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

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

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        // Initial conditions from the original paper (<-- WORKS)
        *((real * )((char *) sv + pitch * 0) + threadID) = -74.7890522727;
        *((real * )((char *) sv + pitch * 1) + threadID) = 136.9896086978;
        *((real * )((char *) sv + pitch * 2) + threadID) = 8.5447311020;
        *((real * )((char *) sv + pitch * 3) + threadID) = 0.0001720623;
        *((real * )((char *) sv + pitch * 4) + threadID) = 0.0184308075;
        *((real * )((char *) sv + pitch * 5) + threadID) = 0.4663168269;
        *((real * )((char *) sv + pitch * 6) + threadID) = 0.3657472179;
        *((real * )((char *) sv + pitch * 7) + threadID) = 0.0486609588;
        *((real * )((char *) sv + pitch * 8) + threadID) = 0.0145766758;
        *((real * )((char *) sv + pitch * 9) + threadID) = 0.2979720207;
        *((real * )((char *) sv + pitch * 10) + threadID) = 0.0692509548;
        *((real * )((char *) sv + pitch * 11) + threadID) = 0.0006146554;
        *((real * )((char *) sv + pitch * 12) + threadID) = 0.0001356656;
        *((real * )((char *) sv + pitch * 13) + threadID) = 0.5943228461;
        *((real * )((char *) sv + pitch * 14) + threadID) = 0.8265709174;
        *((real * )((char *) sv + pitch * 15) + threadID) = 0.9767040566;
        *((real * )((char *) sv + pitch * 16) + threadID) = 0.9717098312;
        *((real * )((char *) sv + pitch * 17) + threadID) = 0.0006830833;
        *((real * )((char *) sv + pitch * 18) + threadID) = 3.2830723338;
        *((real * )((char *) sv + pitch * 19) + threadID) = 0.8199969443;

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

            // Forward Euler variables
            *((real *)((char *)sv) + sv_id) = dt * rDY[0] + *((real *)((char *)sv) + sv_id);
            *((real *)((char *)sv + pitch * 1) + sv_id) = dt * rDY[1] + *((real *)((char *)sv + pitch * 1) + sv_id);
            *((real *)((char *)sv + pitch * 2) + sv_id) = dt * rDY[2] + *((real *)((char *)sv + pitch * 2) + sv_id);
            *((real *)((char *)sv + pitch * 3) + sv_id) = dt * rDY[3] + *((real *)((char *)sv + pitch * 3) + sv_id);
            *((real *)((char *)sv + pitch * 11) + sv_id) = dt * rDY[11] + *((real *)((char *)sv + pitch * 11) + sv_id);
            *((real *)((char *)sv + pitch * 18) + sv_id) = dt * rDY[18] + *((real *)((char *)sv + pitch * 18) + sv_id);
            *((real *)((char *)sv + pitch * 19) + sv_id) = dt * rDY[19] + *((real *)((char *)sv + pitch * 19) + sv_id);

            // Rush Larsen variables
            *((real *)((char *)sv + pitch * 4) + sv_id) = rDY[4];
            *((real *)((char *)sv + pitch * 5) + sv_id) = rDY[5];
            *((real *)((char *)sv + pitch * 6) + sv_id) = rDY[6];
            *((real *)((char *)sv + pitch * 7) + sv_id) = rDY[7];
            *((real *)((char *)sv + pitch * 8) + sv_id) = rDY[8];
            *((real *)((char *)sv + pitch * 9) + sv_id) = rDY[9];
            *((real *)((char *)sv + pitch * 10) + sv_id) = rDY[10];
            *((real *)((char *)sv + pitch * 12) + sv_id) = rDY[12];
            *((real *)((char *)sv + pitch * 13) + sv_id) = rDY[13];
            *((real *)((char *)sv + pitch * 14) + sv_id) = rDY[14];
            *((real *)((char *)sv + pitch * 15) + sv_id) = rDY[15];
            *((real *)((char *)sv + pitch * 16) + sv_id) = rDY[16];
            *((real *)((char *)sv + pitch * 17) + sv_id) = rDY[17];           

        }

    }
}

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt) {

    // States variables
    const real V = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real K_i = *((real*)((char*)sv_ + pitch * 1) + threadID_);
    const real Na_i = *((real*)((char*)sv_ + pitch * 2) + threadID_);
    const real Ca_i = *((real*)((char*)sv_ + pitch * 3) + threadID_);
    const real y = *((real*)((char*)sv_ + pitch * 4) + threadID_);
    const real Xr1 = *((real*)((char*)sv_ + pitch * 5) + threadID_);
    const real Xr2 = *((real*)((char*)sv_ + pitch * 6) + threadID_);
    const real Xs = *((real*)((char*)sv_ + pitch * 7) + threadID_);
    const real m = *((real*)((char*)sv_ + pitch * 8) + threadID_);
    const real h = *((real*)((char*)sv_ + pitch * 9) + threadID_);
    const real j = *((real*)((char*)sv_ + pitch * 10) + threadID_);
    const real Ca_ss = *((real*)((char*)sv_ + pitch * 11) + threadID_);
    const real d = *((real*)((char*)sv_ + pitch * 12) + threadID_);
    const real f = *((real*)((char*)sv_ + pitch * 13) + threadID_);
    const real f2 = *((real*)((char*)sv_ + pitch * 14) + threadID_);
    const real fCass = *((real*)((char*)sv_ + pitch * 15) + threadID_);
    const real s = *((real*)((char*)sv_ + pitch * 16) + threadID_);
    const real r = *((real*)((char*)sv_ + pitch * 17) + threadID_);
    const real Ca_SR = *((real*)((char*)sv_ + pitch * 18) + threadID_);
    const real R_prime = *((real*)((char*)sv_ + pitch * 19) + threadID_);

    // This statement if to avoid instability problems when we have a transmembrane potential below -70mV, 
    // which generates NaN on the solution from the ODEs
    //if (V < INITIAL_V)
    //    V = INITIAL_V;

    // Constants
    const real R = 8314.472;
    const real T = 310;
    const real F = 96485.3415;
    const real Cm = 0.185;
    const real V_c = 0.016404;
    const real P_kna = 0.03;
    const real K_o = 5.4;
    const real Na_o = 140;
    const real Ca_o = 2;
    const real g_f_Na = 0.0145654;
    const real g_f_K = 0.0234346;
    const real g_K1 = 0.065;
    const real g_Kr = 0.0918;
    const real g_Ks = 0.2352;
    const real g_Na = 130.5744;
    const real g_bna = 0.00029;
    const real g_CaL = 3.98e-5;
    const real g_bca = 0.000592;
    const real g_to = 0.08184;
    const real g_sus = 0.0227;
    const real P_NaK = 2.724;
    const real K_mK = 1;
    const real K_mNa = 40;
    const real K_NaCa = 1000;
    const real K_sat = 0.1;
    const real alpha = 2.5;
    const real gamma = 0.35;
    const real Km_Ca = 1.38;
    const real Km_Nai = 87.5;
    const real g_pCa = 0.1238;
    const real K_pCa = 0.0005;
    const real g_pK = 0.0146;
    const real k1_prime = 0.15;
    const real k2_prime = 0.045;
    const real k3 = 0.06;
    const real k4 = 0.005;
    const real EC = 1.5;
    const real max_sr = 2.5;
    const real min_sr = 1;
    const real V_rel = 0.102;
    const real V_xfer = 0.0038;
    const real K_up = 0.00025;
    const real V_leak = 0.00036;
    const real Vmax_up = 0.006375;
    const real Buf_c = 0.2;
    const real K_buf_c = 0.001;
    const real Buf_sr = 10;
    const real K_buf_sr = 0.3;
    const real Buf_ss = 0.4;
    const real K_buf_ss = 0.00025;
    const real V_sr = 0.001094;
    const real V_ss = 5.468e-5;

    // Algebraics
    real f_inf = 1.00000/(1.00000+exp((V+20.0000)/7.00000));
    real tau_f =  1102.50*exp(- powf(V+27.0000, 2.00000)/225.000)+200.000/(1.00000+exp((13.0000 - V)/10.0000))+180.000/(1.00000+exp((V+30.0000)/10.0000))+20.0000;
    real f2_inf = 0.670000/(1.00000+exp((V+35.0000)/7.00000))+0.330000;
    real tau_f2 =  562.000*exp(- powf(V+27.0000, 2.00000)/240.000)+31.0000/(1.00000+exp((25.0000 - V)/10.0000))+80.0000/(1.00000+exp((V+30.0000)/10.0000));
    real fCass_inf = 0.600000/(1.00000+powf(Ca_ss/0.0500000, 2.00000))+0.400000;
    real tau_fCass = 80.0000/(1.00000+powf(Ca_ss/0.0500000, 2.00000))+2.00000;
    real s_inf = 1.00000/(1.00000+exp((V+27.0000)/13.0000));
    real tau_s =  85.0000*exp(- powf(V+25.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((V - 40.0000)/5.00000))+42.0000;
    real r_inf = 1.00000/(1.00000+exp((20.0000 - V)/13.0000));
    real tau_r =  10.4500*exp(- powf(V+40.0000, 2.00000)/1800.00)+7.30000;
    real y_inf = 1.00000/(1.00000+exp((V+80.6000)/6.80000));
    real alpha_y =  1.00000*exp(- 2.90000 -  0.0400000*V);
    real beta_y =  1.00000*exp(3.60000+ 0.110000*V);
    real tau_y = 4000.00/(alpha_y+beta_y);
    real xr1_inf = 1.00000/(1.00000+exp((- 26.0000 - V)/7.00000));
    real alpha_xr1 = 450.000/(1.00000+exp((- 45.0000 - V)/10.0000));
    real beta_xr1 = 6.00000/(1.00000+exp((V+30.0000)/11.5000));
    real tau_xr1 =  1.00000*alpha_xr1*beta_xr1;
    real xr2_inf = 1.00000/(1.00000+exp((V+88.0000)/24.0000));
    real alpha_xr2 = 3.00000/(1.00000+exp((- 60.0000 - V)/20.0000));
    real beta_xr2 = 1.12000/(1.00000+exp((V - 60.0000)/20.0000));
    real tau_xr2 =  1.00000*alpha_xr2*beta_xr2;
    real xs_inf = 1.00000/(1.00000+exp((- 5.00000 - V)/14.0000));
    real alpha_xs = 1400.00/ powf((1.00000+exp((5.00000 - V)/6.00000)), 1.0 / 2);
    real beta_xs = 1.00000/(1.00000+exp((V - 35.0000)/15.0000));
    real tau_xs =  1.00000*alpha_xs*beta_xs+80.0000;
    real m_inf = 1.00000/powf(1.00000+exp((- 56.8600 - V)/9.03000), 2.00000);
    real alpha_m = 1.00000/(1.00000+exp((- 60.0000 - V)/5.00000));
    real beta_m = 0.100000/(1.00000+exp((V+35.0000)/5.00000))+0.100000/(1.00000+exp((V - 50.0000)/200.000));
    real tau_m =  1.00000*alpha_m*beta_m;
    real h_inf = 1.00000/powf(1.00000+exp((V+71.5500)/7.43000), 2.00000);
    real alpha_h = (V<- 40.0000 ?  0.0570000*exp(- (V+80.0000)/6.80000) : 0.00000);
    real beta_h = (V<- 40.0000 ?  2.70000*exp( 0.0790000*V)+ 310000.*exp( 0.348500*V) : 0.770000/( 0.130000*(1.00000+exp((V+10.6600)/- 11.1000))));
    real tau_h = 1.00000/(alpha_h+beta_h);
    real j_inf = 1.00000/powf(1.00000+exp((V+71.5500)/7.43000), 2.00000);
    real alpha_j = (V<- 40.0000 ? (( ( - 25428.0*exp( 0.244400*V) -  6.94800e-06*exp( - 0.0439100*V))*(V+37.7800))/1.00000)/(1.00000+exp( 0.311000*(V+79.2300))) : 0.00000);
    real beta_j = (V<- 40.0000 ? ( 0.0242400*exp( - 0.0105200*V))/(1.00000+exp( - 0.137800*(V+40.1400))) : ( 0.600000*exp( 0.0570000*V))/(1.00000+exp( - 0.100000*(V+32.0000))));
    real tau_j = 1.00000/(alpha_j+beta_j);
    real d_inf = 1.00000/(1.00000+exp((- 8.00000 - V)/7.50000));
    real alpha_d = 1.40000/(1.00000+exp((- 35.0000 - V)/13.0000))+0.250000;
    real beta_d = 1.40000/(1.00000+exp((V+5.00000)/5.00000));
    real gamma_d = 1.00000/(1.00000+exp((50.0000 - V)/20.0000));
    real tau_d =  1.00000*alpha_d*beta_d+gamma_d;
    real i_NaK = (( (( P_NaK*K_o)/(K_o+K_mK))*Na_i)/(Na_i+K_mNa))/(1.00000+ 0.124500*exp(( - 0.100000*V*F)/( R*T))+ 0.0353000*exp(( - V*F)/( R*T)));
    real E_Na =  (( R*T)/F)*log(Na_o/Na_i);
    real i_Na =  g_Na*powf(m, 3.00000)*h*j*(V - E_Na);
    real i_b_Na =  g_bna*(V - E_Na);
    real i_NaCa = ( K_NaCa*( exp(( gamma*V*F)/( R*T))*powf(Na_i, 3.00000)*Ca_o -  exp(( (gamma - 1.00000)*V*F)/( R*T))*powf(Na_o, 3.00000)*Ca_i*alpha))/( (powf(Km_Nai, 3.00000)+powf(Na_o, 3.00000))*(Km_Ca+Ca_o)*(1.00000+ K_sat*exp(( (gamma - 1.00000)*V*F)/( R*T))));
    real i_f_Na =  y*g_f_Na*(V - E_Na);
    real E_K =  (( R*T)/F)*log(K_o/K_i);
    real xK1_inf = 1.00000/(1.00000+exp( 0.100000*(V+75.4400)));
    real i_K1 =  g_K1*xK1_inf*((V - 8.00000) - E_K);
    real i_to =  g_to*r*s*(V - E_K);
    real a = 1.00000/(1.00000+exp((5.00000 - V)/17.0000));
    real i_sus =  g_sus*a*(V - E_K);
    real i_Kr =  g_Kr* powf((K_o/5.40000), 1.0 / 2)*Xr1*Xr2*(V - E_K);
    real E_Ks =  (( R*T)/F)*log((K_o+ P_kna*Na_o)/(K_i+ P_kna*Na_i));
    real i_Ks =  g_Ks*powf(Xs, 2.00000)*(V - E_Ks);
    real i_CaL = ( (( g_CaL*d*f*f2*fCass*4.00000*(V - 15.0000)*powf(F, 2.00000))/( R*T))*( 0.250000*Ca_ss*exp(( 2.00000*(V - 15.0000)*F)/( R*T)) - Ca_o))/(exp(( 2.00000*(V - 15.0000)*F)/( R*T)) - 1.00000);
    real E_Ca =  (( 0.500000*R*T)/F)*log(Ca_o/Ca_i);
    real i_b_Ca =  g_bca*(V - E_Ca);
    real i_p_K = ( g_pK*(V - E_K))/(1.00000+exp((25.0000 - V)/5.98000));
    real i_p_Ca = ( g_pCa*Ca_i)/(Ca_i+K_pCa);
    real i_f_K =  y*g_f_K*(V - E_K);
    real i_f_in = i_f_Na+i_f_K;
    real i_up = Vmax_up/(1.00000+powf(K_up, 2.00000)/powf(Ca_i, 2.00000));
    real i_leak =  V_leak*(Ca_SR - Ca_i);
    real i_xfer =  V_xfer*(Ca_ss - Ca_i);
    real Ca_i_bufc = 1.00000/(1.00000+( Buf_c*K_buf_c)/powf(Ca_i+K_buf_c, 2.00000));
    real kcasr = max_sr - (max_sr - min_sr)/(1.00000+powf(EC/Ca_SR, 2.00000));
    real k2 =  k2_prime*kcasr;
    real k1 = k1_prime/kcasr;
    real O = ( k1*powf(Ca_ss, 2.00000)*R_prime)/(k3+ k1*powf(Ca_ss, 2.00000));
    real i_rel =  V_rel*O*(Ca_SR - Ca_ss);
    real Ca_sr_bufsr = 1.00000/(1.00000+( Buf_sr*K_buf_sr)/powf(Ca_SR+K_buf_sr, 2.00000));
    real Ca_ss_bufss = 1.00000/(1.00000+( Buf_ss*K_buf_ss)/powf(Ca_ss+K_buf_ss, 2.00000));

    //  ** I manually added the stimulus current
    // Rates
    real d_dt_V = (- 1.00000/1.00000)*(i_K1+i_to+i_sus+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_b_Na+i_NaCa+i_b_Ca+i_p_K+i_p_Ca+i_f_in+stim_current);

    real d_dt_K_i = (( - 1.00000*((i_K1+i_to+i_f_K+i_sus+i_Kr+i_Ks+i_p_K) -  2.00000*i_NaK))/( 1.00000*V_c*F))*Cm;
    real d_dt_Na_i = (( - 1.00000*(i_Na+i_b_Na+i_f_Na+ 3.00000*i_NaK+ 3.00000*i_NaCa))/( 1.00000*V_c*F))*Cm;
    real d_dt_Ca_i = Ca_i_bufc*((( (i_leak - i_up)*V_sr)/V_c+i_xfer) - ( 1.00000*((i_b_Ca+i_p_Ca) -  2.00000*i_NaCa)*Cm)/( 2.00000*1.00000*V_c*F));
    real d_dt_Ca_ss = Ca_ss_bufss*((( - 1.00000*i_CaL*Cm)/( 2.00000*1.00000*V_ss*F)+( i_rel*V_sr)/V_ss) - ( i_xfer*V_c)/V_ss);
    real d_dt_Ca_SR = Ca_sr_bufsr*(i_up - (i_rel+i_leak));
    real d_dt_R_prime = - k2*Ca_ss*R_prime+ k4*(1.00000 - R_prime);

    // Forward Euler variables
    rDY_[0] = d_dt_V;
    rDY_[1] = d_dt_K_i;
    rDY_[2] = d_dt_Na_i;
    rDY_[3] = d_dt_Ca_i;
    rDY_[11] = d_dt_Ca_ss;
    rDY_[18] = d_dt_Ca_SR;
    rDY_[19] = d_dt_R_prime; 

    // Rush Larsen variables
    rDY_[4] = y_inf + (y-y_inf)*expf(-dt/tau_y);
    rDY_[5] = xr1_inf + (Xr1-xr1_inf)*expf(-dt/tau_xr1);
    rDY_[6] = xr2_inf + (Xr2-xr2_inf)*expf(-dt/tau_xr2);
    rDY_[7] = xs_inf + (Xs-xs_inf)*expf(-dt/tau_xs);
    rDY_[8] = m_inf + (m-m_inf)*expf(-dt/tau_m);
    rDY_[9] = h_inf + (h-h_inf)*expf(-dt/tau_h);
    rDY_[10] = j_inf + (j-j_inf)*expf(-dt/tau_j);
    rDY_[12] = d_inf + (d-d_inf)*expf(-dt/tau_d);
    rDY_[13] = f_inf + (f-f_inf)*expf(-dt/tau_f);
    rDY_[14] = f2_inf + (f2-f2_inf)*expf(-dt/tau_f2);
    rDY_[15] = fCass_inf + (fCass-fCass_inf)*expf(-dt/tau_fCass);
    rDY_[16] = s_inf + (s-s_inf)*expf(-dt/tau_s);
    rDY_[17] = r_inf + (r-r_inf)*expf(-dt/tau_r);

}

