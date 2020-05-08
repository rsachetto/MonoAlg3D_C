#include "noble_1962.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    // Default CellML 
    //sv[0] = -75.5344986658f;    // V millivolt
    //sv[1] = 0.060546727200f;    // m dimensionless
    //sv[2] = 0.725900135500f;    // h dimensionless
    //sv[3] = 0.470923970800f;    // n dimensionless

    // BCL = 300ms | 10 pulses
    sv[0] = -81.1893;    // V millivolt
    sv[1] = 0.0443563;    // m dimensionless
    sv[2] = 0.851652;    // h dimensionless
    sv[3] = 0.58291;    // n dimensionless

    // BCL = 500ms | 30 pulses
    //sv[0] = -75.238;
    //sv[1] = 0.0615111;
    //sv[2] = 0.718401;
    //sv[3] = 0.467409;

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

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt);

    //THIS MODEL USES THE Rush Larsen Method TO SOLVE THE EDOS
    // Forward Euler variables
    sv[0] = dt*rDY[0] + rY[0];

    // Rush Larsen variables
    sv[1] = rDY[1];
    sv[2] = rDY[2];
    sv[3] = rDY[3];
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) 
{

    //State variables
    const real V_old_ = sv[0];
    const real m_old_ = sv[1];
    const real h_old_ = sv[2];
    const real n_old_ = sv[3];

    //Parameters
    //const real Cm = 12.00000000000000000e+00f;             // (microF)
    //const real g_na_max = 400000.00000000000000000e+00f;   // (microS)
    //const real E_na = 40.00000000000000000e+00f;           // (millivolt)
    //const real g_L = 75.00000000000000000e+00f;            // (microS)
    //const real E_L = -60.00000000000000000e+00f;           // (millivolt)
    const real Cm = 12.0;             // (microF)
    const real g_na_max = 400.0;   // (microS)
    const real E_na = 40.0;           // (millivolt)
    const real g_L = 0.075;            // (microS)
    const real E_L = -60.0;           // (millivolt)


    real calc_I_stim = stim_current;

    // Algebraics
    //real g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    //real alpha_m = ( 100.000*(- V_old_ - 48.0000))/(exp((- V_old_ - 48.0000)/15.0000) - 1.00000);
    //real alpha_h =  170.000*exp((- V_old_ - 90.0000)/20.0000);
    //real alpha_n = ( 0.100000*(- V_old_ - 50.0000))/(exp((- V_old_ - 50.0000)/10.0000) - 1.00000);
    //real i_na =  (g_na+140.000)*(V_old_ - E_na);
    //real i_na_no_oscilation = (g_na+122.500)*(V_old_ - E_na);
    //real beta_m = ( 120.000*(V_old_+8.00000))/(exp((V_old_+8.00000)/5.00000) - 1.00000);
    //real beta_h = 1000.00/(1.00000+exp((- V_old_ - 42.0000)/10.0000));
    //real beta_n =  2.00000*exp((- V_old_ - 90.0000)/80.0000);
    //real g_K1 =  1200.00*exp((- V_old_ - 90.0000)/50.0000)+ 15.0000*exp((V_old_+90.0000)/60.0000);
    //real g_K2 =  1200.00*pow(n_old_, 4.00000);
    //real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    //real i_leak =  g_L*(V_old_ - E_L);

    real g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    real alpha_h = ((1.7e-01*exp((((-V_old_)-9.0e+01)/2.0e+01))));
    real alpha_m = (((1.0e-01*((-V_old_)-4.8e+01))/(exp((((-V_old_)-4.8e+01)/1.5e+01))-1.0e+00)));
    real alpha_n = (((1.0e-04*((-V_old_)-5.0e+01))/(exp((((-V_old_)-5.0e+01)/1.0e+01))-1.0e+00)));
    real i_na = (g_na+1.4e-01)*(V_old_ - E_na);
    real i_na_no_oscilation = (g_na+1.2e-01)*(V_old_ - E_na);
    real beta_m = (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
    real beta_h = ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
    real beta_n = ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
    real g_K1 = 1.2*exp((((-V_old_)-9.0e+01)/5.0e+01)) + (1.5e-02*exp(((V_old_+9.0e+01)/6.0e+01)));
    real g_K2 = 1.2*pow(n_old_,4.0e+00);
    real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    real i_leak =  g_L*(V_old_ - E_L);
    real tau_h = 1.0 / (alpha_h + beta_h);
    real tau_m = 1.0 / (alpha_m + beta_m);
    real tau_n = 1.0 / (alpha_n + beta_n);
    real inf_h = alpha_h / (alpha_h + beta_h);
    real inf_m = alpha_m / (alpha_m + beta_m);
    real inf_n = alpha_n / (alpha_n + beta_n);

    // Rates
    //rDY_[0] = (- (i_na + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    //rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  beta_m*m_old_) * 1.0E-03;
    //rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  beta_h*h_old_) * 1.0E-03;
    //rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  beta_n*n_old_) * 1.0E-03;

    rDY_[0] = (- (i_na + i_k + i_leak + calc_I_stim)/Cm);
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm);
    //rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  beta_m*m_old_);
    //rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  beta_h*h_old_);
    //rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  beta_n*n_old_);

    // Rush Larsen
    rDY_[1] = inf_m + (m_old_-inf_m)*expf(-dt/tau_m);
    rDY_[2] = inf_h + (h_old_-inf_h)*expf(-dt/tau_h);
    rDY_[3] = inf_n + (n_old_-inf_n)*expf(-dt/tau_n);

}
