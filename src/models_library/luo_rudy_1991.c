#include "luo_rudy_1991.h"

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    sv[0] = -84.380111f; //V millivolt 
    sv[1] = 0.001713f; //m dimensionless 
    sv[2] = 0.982661f; //h dimensionless 
    sv[3] = 0.989108f; //j dimensionless 
    sv[4] = 0.003021f; //d dimensionless 
    sv[5] = 0.999968f; //f dimensionless 
    sv[6] = 0.041760f; //X dimensionless 
    sv[7] = 0.000179f; //Cai millimolar 
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

    RHS_cpu(rY, rDY, stim_current);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

#define IFNUMBER_1(name)if((V_old_<(-4.000000000000000e+01f))) { (name) = (1.350000000000000e-01f*expf(((8.000000000000000e+01f+V_old_)/(-6.800000000000000e+00f))));    }  else{ (name) = 0.000000000000000e+00f;    }
#define IFNUMBER_2(name)if((V_old_<(-4.000000000000000e+01f))) { (name) = ((3.560000000000000e+00f*expf((7.900000000000000e-02f*V_old_)))+(3.100000000000000e+05f*expf((3.500000000000000e-01f*V_old_))));    }  else{ (name) = (1.000000000000000e+00f/(1.300000000000000e-01f*(1.000000000000000e+00f+expf(((V_old_+1.066000000000000e+01f)/(-1.110000000000000e+01f))))));    }
#define IFNUMBER_3(name)if((V_old_<(-4.000000000000000e+01f))) { (name) = (((((-1.271400000000000e+05f)*expf((2.444000000000000e-01f*V_old_)))-(3.474000000000000e-05f*expf(((-4.391000000000000e-02f)*V_old_))))*(V_old_+3.778000000000000e+01f))/(1.000000000000000e+00f+expf((3.110000000000000e-01f*(V_old_+7.923000000000000e+01f)))));    }  else{ (name) = 0.000000000000000e+00f;    }
#define IFNUMBER_4(name)if((V_old_<(-4.000000000000000e+01f))) { (name) = ((1.212000000000000e-01f*expf(((-1.052000000000000e-02f)*V_old_)))/(1.000000000000000e+00f+expf(((-1.378000000000000e-01f)*(V_old_+4.014000000000000e+01f)))));    }  else{ (name) = ((3.000000000000000e-01f*expf(((-2.535000000000000e-07f)*V_old_)))/(1.000000000000000e+00f+expf(((-1.000000000000000e-01f)*(V_old_+3.200000000000000e+01f)))));    }
#define IFNUMBER_5(name)if((V_old_>(-1.000000000000000e+02f))) { (name) = ((2.837000000000000e+00f*(expf((4.000000000000000e-02f*(V_old_+7.700000000000000e+01f)))-1.000000000000000e+00f))/((V_old_+7.700000000000000e+01f)*expf((4.000000000000000e-02f*(V_old_+3.500000000000000e+01f)))));    }  else{ (name) = 1.000000000000000e+00f;    }

void RHS_cpu(const real *sv, real *rDY_, real stim_current) {

    //State variables
    const real V_old_ = sv[0];
    const real m_old_ = sv[1];
    const real h_old_ = sv[2];
    const real j_old_ = sv[3];
    const real d_old_ = sv[4];
    const real f_old_ = sv[5];
    const real X_old_ = sv[6];
    const real Cai_old_ = sv[7];

    //Parameters
    const real C = 1.000000000000000e+00f;
    const real R = 8.314000000000000e+03f;
    const real T = 3.100000000000000e+02f;
    const real F = 9.648460000000001e+04f;
    const real Nao = 1.400000000000000e+02f;
    const real Nai = 1.800000000000000e+01f;
    const real g_Na = 2.300000000000000e+01f;
    const real Ko = 5.400000000000000e+00f;
    const real PR_NaK = 1.833000000000000e-02f;
    const real Ki = 1.450000000000000e+02f;
    const real g_Kp = 1.830000000000000e-02f;
    const real g_b = 3.921000000000000e-02f;
    const real E_b = -5.987000000000000e+01f;

    real calc_I_stim = stim_current;
    real calc_E_Na = (((R*T)/F)*logf((Nao/Nai)));	//2
    real calc_alpha_m = ((3.200000000000000e-01f*(V_old_+4.713000000000000e+01f))/(1.000000000000000e+00f-expf(((-1.000000000000000e-01f)*(V_old_+4.713000000000000e+01f)))));	//4
    real calc_beta_m = (8.000000000000000e-02f*expf(((-V_old_)/1.100000000000000e+01f)));	//5
    real calc_alpha_h = 0.0f;
    IFNUMBER_1(calc_alpha_h);	//7
    real calc_beta_h = 0.0f;
    IFNUMBER_2(calc_beta_h);	//8
    real calc_alpha_j = 0.0f;
    IFNUMBER_3(calc_alpha_j);	//10
    real calc_beta_j = 0.0f;
    IFNUMBER_4(calc_beta_j);	//11
    real calc_E_si = (7.700000000000000e+00f-(1.302870000000000e+01f*logf((Cai_old_/1.000000000000000e+00f))));	//13
    real calc_alpha_d = ((9.500000000000000e-02f*expf(((-1.000000000000000e-02f)*(V_old_-5.000000000000000e+00f))))/(1.000000000000000e+00f+expf(((-7.199999999999999e-02f)*(V_old_-5.000000000000000e+00f)))));	//15
    real calc_beta_d = ((7.000000000000001e-02f*expf(((-1.700000000000000e-02f)*(V_old_+4.400000000000000e+01f))))/(1.000000000000000e+00f+expf((5.000000000000000e-02f*(V_old_+4.400000000000000e+01f)))));	//16
    real calc_alpha_f = ((1.200000000000000e-02f*expf(((-8.000000000000000e-03f)*(V_old_+2.800000000000000e+01f))))/(1.000000000000000e+00f+expf((1.500000000000000e-01f*(V_old_+2.800000000000000e+01f)))));	//18
    real calc_beta_f = ((6.500000000000000e-03f*expf(((-2.000000000000000e-02f)*(V_old_+3.000000000000000e+01f))))/(1.000000000000000e+00f+expf(((-2.000000000000000e-01f)*(V_old_+3.000000000000000e+01f)))));	//19
    real calc_g_K = (2.820000000000000e-01f*powf((Ko/5.400000000000000e+00f),0.5f));	//21
    real calc_E_K = (((R*T)/F)*logf(((Ko+(PR_NaK*Nao))/(Ki+(PR_NaK*Nai)))));	//22
    real calc_alpha_X = ((5.000000000000000e-04f*expf((8.300000000000000e-02f*(V_old_+5.000000000000000e+01f))))/(1.000000000000000e+00f+expf((5.700000000000000e-02f*(V_old_+5.000000000000000e+01f)))));	//24
    real calc_beta_X = ((1.300000000000000e-03f*expf(((-6.000000000000000e-02f)*(V_old_+2.000000000000000e+01f))))/(1.000000000000000e+00f+expf(((-4.000000000000000e-02f)*(V_old_+2.000000000000000e+01f)))));	//25
    real calc_Xi = 0.0f;
    IFNUMBER_5(calc_Xi);	//27
    real calc_g_K1 = (6.047000000000000e-01f*powf((Ko/5.400000000000000e+00f),0.5f));	//28
    real calc_E_K1 = (((R*T)/F)*logf((Ko/Ki)));	//29
    real calc_Kp = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((7.488000000000000e+00f-V_old_)/5.980000000000000e+00f))));	//35
    real calc_i_b = (g_b*(V_old_-E_b));	//37
    real calc_i_Na = (g_Na*powf(m_old_,3.000000000000000e+00f)*h_old_*j_old_*(V_old_-calc_E_Na));	//3
    real calc_i_si = (9.000000000000000e-02f*d_old_*f_old_*(V_old_-calc_E_si));	//14
    real calc_alpha_K1 = (1.020000000000000e+00f/(1.000000000000000e+00f+expf((2.385000000000000e-01f*((V_old_-calc_E_K1)-5.921500000000000e+01f)))));	//31
    real calc_beta_K1 = (((4.912400000000000e-01f*expf((8.032000000000000e-02f*((V_old_+5.476000000000000e+00f)-calc_E_K1))))+(1.000000000000000e+00f*expf((6.175000000000000e-02f*(V_old_-(calc_E_K1+5.943099999999999e+02f))))))/(1.000000000000000e+00f+expf(((-5.143000000000000e-01f)*((V_old_-calc_E_K1)+4.753000000000000e+00f)))));	//32
    real calc_E_Kp = calc_E_K1;	//34
    real calc_i_K = (calc_g_K*X_old_*calc_Xi*(V_old_-calc_E_K));	//23
    real calc_K1_infinity = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));	//33
    real calc_i_Kp = (g_Kp*calc_Kp*(V_old_-calc_E_Kp));	//36
    real calc_i_K1 = (calc_g_K1*calc_K1_infinity*(V_old_-calc_E_K1));	//30

    rDY_[0] = (((-1.000000000000000e+00f)/C)*(calc_I_stim+calc_i_Na+calc_i_si+calc_i_K+calc_i_K1+calc_i_Kp+calc_i_b));
    rDY_[1] = ((calc_alpha_m*(1.000000000000000e+00f-m_old_))-(calc_beta_m*m_old_));
    rDY_[2] = ((calc_alpha_h*(1.000000000000000e+00f-h_old_))-(calc_beta_h*h_old_));
    rDY_[3] = ((calc_alpha_j*(1.000000000000000e+00f-j_old_))-(calc_beta_j*j_old_));
    rDY_[4] = ((calc_alpha_d*(1.000000000000000e+00f-d_old_))-(calc_beta_d*d_old_));
    rDY_[5] = ((calc_alpha_f*(1.000000000000000e+00f-f_old_))-(calc_beta_f*f_old_));
    rDY_[6] = ((calc_alpha_X*(1.000000000000000e+00f-X_old_))-(calc_beta_X*X_old_));
    rDY_[7] = ((((-1.000000000000000e-04f)/1.000000000000000e+00f)*calc_i_si)+(7.000000000000001e-02f*(1.000000000000000e-04f-Cai_old_)));

}

