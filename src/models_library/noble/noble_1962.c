#include "noble_1962.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    // Normal
    //sv[0] = -75.5344986658f;    // V millivolt
    //sv[1] = 0.060546727200f;    // m dimensionless
    //sv[2] = 0.725900135500f;    // h dimensionless
    //sv[3] = 0.470923970800f;    // n dimensionless

    uint32_t num_volumes = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_volumes*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_volumes; i++) {
        real *sv = &solver->sv[i * NEQ];
        // BCL = 300ms
        sv[0] = -81.1893;  // V millivolt
        sv[1] = 0.0443563; // m dimensionless
        sv[2] = 0.851652;  // h dimensionless
        sv[3] = 0.58291;   // n dimensionless
    }

}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    uint32_t sv_id;

    OMP(parallel for private(sv_id))
    for (uint32_t i = 0; i < num_cells_to_solve; i++) {

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

    RHS_cpu(rY, rDY, stim_current);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current) {

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
    real alpha_m = (((1.0e-01*((-V_old_)-4.8e+01))/(expm1((((-V_old_)-4.8e+01)/1.5e+01)))));
    real alpha_n = (((1.0e-04*((-V_old_)-5.0e+01))/(expm1((((-V_old_)-5.0e+01)/1.0e+01)))));
    real i_na = (g_na+1.4e-01)*(V_old_ - E_na);
//    real i_na_no_oscilation = (g_na+1.2e-01)*(V_old_ - E_na);
    double beta_m = (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
    double beta_h = ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
    double beta_n = ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
    real g_K1 = 1.2*exp((((-V_old_)-9.0e+01)/5.0e+01)) + (1.5e-02*exp(((V_old_+9.0e+01)/6.0e+01)));
    real g_K2 = 1.2*pow(n_old_,4.0e+00);
    real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    real i_leak =  g_L*(V_old_ - E_L);

    // Rates
    //rDY_[0] = (- (i_na + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    //rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  beta_m*m_old_) * 1.0E-03;
    //rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  beta_h*h_old_) * 1.0E-03;
    //rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  beta_n*n_old_) * 1.0E-03;

    rDY_[0] = (- (i_na + i_k + i_leak + calc_I_stim)/Cm);
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm);
    rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  beta_m*m_old_);
    rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  beta_h*h_old_);
    rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  beta_n*n_old_);

}

