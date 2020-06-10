#include <stdio.h>
#include "stewart_aslanidi_noble_2009.h"

GET_CELL_MODEL_DATA(init_cell_model_data) 
{

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) 
{

    log_to_stdout_and_file("Using Stewart-Aslanidi-Noble 2009 CPU model\n");

    uint32_t num_volumes = solver->original_num_cells;
	solver->sv = (real*)malloc(NEQ*num_volumes*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_volumes; i++) {
        real *sv = &solver->sv[i * NEQ];

        // Initial conditions from the original paper (<-- WORKS)
        sv[0] = -74.7890522727;
        sv[1] = 136.9896086978;
        sv[2] = 8.5447311020;
        sv[3] = 0.0001720623;
        sv[4] = 0.0184308075;
        sv[5] = 0.4663168269;
        sv[6] = 0.3657472179;
        sv[7] = 0.0486609588;
        sv[8] = 0.0145766758;
        sv[9] = 0.2979720207;
        sv[10] = 0.0692509548;
        sv[11] = 0.0006146554;
        sv[12] = 0.0001356656;
        sv[13] = 0.5943228461;
        sv[14] = 0.8265709174;
        sv[15] = 0.9767040566;
        sv[16] = 0.9717098312;
        sv[17] = 0.0006830833;
        sv[18] = 3.2830723338;
        sv[19] = 0.8199969443;
    }

}

SOLVE_MODEL_ODES(solve_model_odes_cpu)
{

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

//    uint32_t *mapping = ((uint32_t*)extra_data);

    OMP(parallel for private(sv_id))
    for (uint32_t i = 0; i < num_cells_to_solve; i++)
    {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = (uint32_t )i;

        for (int j = 0; j < num_steps; ++j) 
        {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);
        }
    }

}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  
{
    real rY[NEQ], rDY[NEQ];

    // Save old value of the state vector
    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    // Solve Right-hand-side of the ODE's
    RHS_cpu(rY, rDY, stim_current,dt);

    // Forward Euler variables
    sv[0] = dt*rDY[0] + rY[0];

    sv[1] = dt*rDY[1] + rY[1];
    sv[2] = dt*rDY[2] + rY[2];
    sv[3] = dt*rDY[3] + rY[3];
    sv[11] = dt*rDY[11] + rY[11];
    sv[18] = dt*rDY[18] + rY[18];
    sv[19] = dt*rDY[19] + rY[19];

    // Rush Larsen variables
    sv[4] = rDY[4];
    sv[5] = rDY[5];
    sv[6] = rDY[6];
    sv[7] = rDY[7];
    sv[8] = rDY[8];
    sv[9] = rDY[9];
    sv[10] = rDY[10];
    sv[12] = rDY[12];
    sv[13] = rDY[13];
    sv[14] = rDY[14];
    sv[15] = rDY[15];
    sv[16] = rDY[16];
    sv[17] = rDY[17];
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) 
{
    // States variables
    const real V = sv[0];
    const real K_i = sv[1];
    const real Na_i = sv[2];
    const real Ca_i = sv[3];
    const real y = sv[4];
    const real Xr1 = sv[5];
    const real Xr2 = sv[6];
    const real Xs = sv[7];
    const real m = sv[8];
    const real h = sv[9];
    const real j = sv[10];
    const real Ca_ss = sv[11];
    const real d = sv[12];
    const real f = sv[13];
    const real f2 = sv[14];
    const real fCass = sv[15];
    const real s = sv[16];
    const real r = sv[17];
    const real Ca_SR = sv[18];
    const real R_prime = sv[19];
    
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
    real tau_f =  1102.50*exp(- pow(V+27.0000, 2.00000)/225.000)+200.000/(1.00000+exp((13.0000 - V)/10.0000))+180.000/(1.00000+exp((V+30.0000)/10.0000))+20.0000;
    real f2_inf = 0.670000/(1.00000+exp((V+35.0000)/7.00000))+0.330000;
    real tau_f2 =  562.000*exp(- pow(V+27.0000, 2.00000)/240.000)+31.0000/(1.00000+exp((25.0000 - V)/10.0000))+80.0000/(1.00000+exp((V+30.0000)/10.0000));
    real fCass_inf = 0.600000/(1.00000+pow(Ca_ss/0.0500000, 2.00000))+0.400000;
    real tau_fCass = 80.0000/(1.00000+pow(Ca_ss/0.0500000, 2.00000))+2.00000;
    real s_inf = 1.00000/(1.00000+exp((V+27.0000)/13.0000));
    real tau_s =  85.0000*exp(- pow(V+25.0000, 2.00000)/320.000)+5.00000/(1.00000+exp((V - 40.0000)/5.00000))+42.0000;
    real r_inf = 1.00000/(1.00000+exp((20.0000 - V)/13.0000));
    real tau_r =  10.4500*exp(- pow(V+40.0000, 2.00000)/1800.00)+7.30000;
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
    real alpha_xs = 1400.00/ pow((1.00000+exp((5.00000 - V)/6.00000)), 1.0 / 2);
    real beta_xs = 1.00000/(1.00000+exp((V - 35.0000)/15.0000));
    real tau_xs =  1.00000*alpha_xs*beta_xs+80.0000;
    real m_inf = 1.00000/pow(1.00000+exp((- 56.8600 - V)/9.03000), 2.00000);
    real alpha_m = 1.00000/(1.00000+exp((- 60.0000 - V)/5.00000));
    real beta_m = 0.100000/(1.00000+exp((V+35.0000)/5.00000))+0.100000/(1.00000+exp((V - 50.0000)/200.000));
    real tau_m =  1.00000*alpha_m*beta_m;
    real h_inf = 1.00000/pow(1.00000+exp((V+71.5500)/7.43000), 2.00000);
    real alpha_h = (V<- 40.0000 ?  0.0570000*exp(- (V+80.0000)/6.80000) : 0.00000);
    real beta_h = (V<- 40.0000 ?  2.70000*exp( 0.0790000*V)+ 310000.*exp( 0.348500*V) : 0.770000/( 0.130000*(1.00000+exp((V+10.6600)/- 11.1000))));
    real tau_h = 1.00000/(alpha_h+beta_h);
    real j_inf = 1.00000/pow(1.00000+exp((V+71.5500)/7.43000), 2.00000);
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
    real i_Na =  g_Na*pow(m, 3.00000)*h*j*(V - E_Na);
    real i_b_Na =  g_bna*(V - E_Na);
    real i_NaCa = ( K_NaCa*( exp(( gamma*V*F)/( R*T))*pow(Na_i, 3.00000)*Ca_o -  exp(( (gamma - 1.00000)*V*F)/( R*T))*pow(Na_o, 3.00000)*Ca_i*alpha))/( (pow(Km_Nai, 3.00000)+pow(Na_o, 3.00000))*(Km_Ca+Ca_o)*(1.00000+ K_sat*exp(( (gamma - 1.00000)*V*F)/( R*T))));
    real i_f_Na =  y*g_f_Na*(V - E_Na);
    real E_K =  (( R*T)/F)*log(K_o/K_i);
    real xK1_inf = 1.00000/(1.00000+exp( 0.100000*(V+75.4400)));
    real i_K1 =  g_K1*xK1_inf*((V - 8.00000) - E_K);
    real i_to =  g_to*r*s*(V - E_K);
    real a = 1.00000/(1.00000+exp((5.00000 - V)/17.0000));
    real i_sus =  g_sus*a*(V - E_K);
    real i_Kr =  g_Kr* pow((K_o/5.40000), 1.0 / 2)*Xr1*Xr2*(V - E_K);
    real E_Ks =  (( R*T)/F)*log((K_o+ P_kna*Na_o)/(K_i+ P_kna*Na_i));
    real i_Ks =  g_Ks*pow(Xs, 2.00000)*(V - E_Ks);
    real i_CaL = ( (( g_CaL*d*f*f2*fCass*4.00000*(V - 15.0000)*pow(F, 2.00000))/( R*T))*( 0.250000*Ca_ss*exp(( 2.00000*(V - 15.0000)*F)/( R*T)) - Ca_o))/(exp(( 2.00000*(V - 15.0000)*F)/( R*T)) - 1.00000);
    real E_Ca =  (( 0.500000*R*T)/F)*log(Ca_o/Ca_i);
    real i_b_Ca =  g_bca*(V - E_Ca);
    real i_p_K = ( g_pK*(V - E_K))/(1.00000+exp((25.0000 - V)/5.98000));
    real i_p_Ca = ( g_pCa*Ca_i)/(Ca_i+K_pCa);
    real i_f_K =  y*g_f_K*(V - E_K);
    real i_f_in = i_f_Na+i_f_K;
    real i_up = Vmax_up/(1.00000+pow(K_up, 2.00000)/pow(Ca_i, 2.00000));
    real i_leak =  V_leak*(Ca_SR - Ca_i);
    real i_xfer =  V_xfer*(Ca_ss - Ca_i);
    real Ca_i_bufc = 1.00000/(1.00000+( Buf_c*K_buf_c)/pow(Ca_i+K_buf_c, 2.00000));
    real kcasr = max_sr - (max_sr - min_sr)/(1.00000+pow(EC/Ca_SR, 2.00000));
    real k2 =  k2_prime*kcasr;
    real k1 = k1_prime/kcasr;
    real O = ( k1*pow(Ca_ss, 2.00000)*R_prime)/(k3+ k1*pow(Ca_ss, 2.00000));
    real i_rel =  V_rel*O*(Ca_SR - Ca_ss);
    real Ca_sr_bufsr = 1.00000/(1.00000+( Buf_sr*K_buf_sr)/pow(Ca_SR+K_buf_sr, 2.00000));
    real Ca_ss_bufss = 1.00000/(1.00000+( Buf_ss*K_buf_ss)/pow(Ca_ss+K_buf_ss, 2.00000));

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

// The automatic pacing from the Purkinje cells can be interrupted by blocking the INa current by 100% (ALGEBRAIC[54] = INa)
