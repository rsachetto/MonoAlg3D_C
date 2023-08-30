
#define IFNUMBER_1(name)if((V_old_==(-4.713000000000000e+01))) { (name) = 3.200000000000000e+00;    }  else{ (name) = ((3.200000000000000e-01*(V_old_+4.713000000000000e+01))/(1.000000000000000e+00-exp(((-1.000000000000000e-01)*(V_old_+4.713000000000000e+01)))));    }
#define IFNUMBER_2(name)if((V_old_<(-4.000000000000000e+01))) { (name) = (1.350000000000000e-01*exp(((V_old_+8.000000000000000e+01)/(-6.800000000000000e+00))));    }  else{ (name) = 0.000000000000000e+00;    }
#define IFNUMBER_3(name)if((V_old_<(-4.000000000000000e+01))) { (name) = ((3.560000000000000e+00*exp((7.900000000000000e-02*V_old_)))+(3.100000000000000e+05*exp((3.500000000000000e-01*V_old_))));    }  else{ (name) = (1.000000000000000e+00/(1.300000000000000e-01*(1.000000000000000e+00+exp(((V_old_+1.066000000000000e+01)/(-1.110000000000000e+01))))));    }
#define IFNUMBER_4(name)if((V_old_<(-4.000000000000000e+01))) { (name) = (((((-1.271400000000000e+05)*exp((2.444000000000000e-01*V_old_)))-(3.474000000000000e-05*exp(((-4.391000000000000e-02)*V_old_))))*(V_old_+3.778000000000000e+01))/(1.000000000000000e+00+exp((3.110000000000000e-01*(V_old_+7.923000000000000e+01)))));    }  else{ (name) = 0.000000000000000e+00;    }
#define IFNUMBER_5(name)if((V_old_<(-4.000000000000000e+01))) { (name) = ((1.212000000000000e-01*exp(((-1.052000000000000e-02)*V_old_)))/(1.000000000000000e+00+exp(((-1.378000000000000e-01)*(V_old_+4.014000000000000e+01)))));    }  else{ (name) = ((3.000000000000000e-01*exp(((-2.535000000000000e-07)*V_old_)))/(1.000000000000000e+00+exp(((-1.000000000000000e-01)*(V_old_+3.200000000000000e+01)))));    }
#define IFNUMBER_6(name)if((fabs((V_old_+1.410000000000000e+01))<1.000000000000000e-10)) { (name) = 1.500000000000000e-03;    }  else{ (name) = ((3.000000000000000e-04*(V_old_+1.410000000000000e+01))/(1.000000000000000e+00-exp(((V_old_+1.410000000000000e+01)/(-5.000000000000000e+00)))));    }
#define IFNUMBER_7(name)if((fabs((V_old_-3.332800000000000e+00))<1.000000000000000e-10)) { (name) = 3.783611800000000e-04;    }  else{ (name) = ((7.389800000000000e-05*(V_old_-3.332800000000000e+00))/(exp(((V_old_-3.332800000000000e+00)/5.123700000000000e+00))-1.000000000000000e+00));    }
#define IFNUMBER_8(name)if((fabs((V_old_-1.990000000000000e+01))<1.000000000000000e-10)) { (name) = 6.800000000000000e-04;    }  else{ (name) = ((4.000000000000000e-05*(V_old_-1.990000000000000e+01))/(1.000000000000000e+00-exp(((V_old_-1.990000000000000e+01)/(-1.700000000000000e+01)))));    }
#define IFNUMBER_9(name)if((fabs((V_old_-1.990000000000000e+01))<1.000000000000000e-10)) { (name) = 3.150000000000000e-04;    }  else{ (name) = ((3.500000000000000e-05*(V_old_-1.990000000000000e+01))/(exp(((V_old_-1.990000000000000e+01)/9.000000000000000e+00))-1.000000000000000e+00));    }
#define IFNUMBER_10(name)if((fabs((V_old_+1.000000000000000e+01))<1.000000000000000e-10)) { (name) = (4.579000000000000e+00/(1.000000000000000e+00+exp(((V_old_+1.000000000000000e+01)/(-6.240000000000000e+00)))));    }  else{ (name) = ((1.000000000000000e+00-exp(((V_old_+1.000000000000000e+01)/(-6.240000000000000e+00))))/(3.500000000000000e-02*(V_old_+1.000000000000000e+01)*(1.000000000000000e+00+exp(((V_old_+1.000000000000000e+01)/(-6.240000000000000e+00))))));    }
#define IFNUMBER_11(name)if((fabs((V_old_-7.900000000000000e+00))<1.000000000000000e-10)) { (name) = ((6.000000000000000e+00*2.000000000000000e-01)/1.300000000000000e+00);    }  else{ (name) = ((6.000000000000000e+00*(1.000000000000000e+00-exp(((-(V_old_-7.900000000000000e+00))/5.000000000000000e+00))))/((1.000000000000000e+00+(3.000000000000000e-01*exp(((-(V_old_-7.900000000000000e+00))/5.000000000000000e+00))))*1.000000000000000e+00*(V_old_-7.900000000000000e+00)));    }    //Parameters

const real Cm = 1.000000000000000e+02f;
const real g_Na = 7.800000000000000e+00f;
const real R = 8.314299999999999e+00f;
const real T = 3.100000000000000e+02f;
const real F = 9.648670000000000e+01f;
const real Na_o = 1.400000000000000e+02f;
const real K_o = 5.400000000000000e+00f;
const real g_K1 = 9.000000000000000e-02f;
const real g_to = 1.652000000000000e-01f;
const real K_Q10 = 3.000000000000000e+00f;
const real g_Kr = 2.941176500000000e-02f;
const real g_Ks = 1.294117600000000e-01f;
const real g_Ca_L = 1.237500000000000e-01f;
const real i_NaK_max = 5.993387400000000e-01f;
const real Km_Na_i = 1.000000000000000e+01f;
const real Km_K_o = 1.500000000000000e+00f;
const real Ca_o = 1.800000000000000e+00f;
const real g_B_Na = 6.744375000000000e-04f;
const real g_B_Ca = 1.131000000000000e-03f;
const real g_B_K = 0.000000000000000e+00f;
const real I_NaCa_max = 1.600000000000000e+03f;
const real gamma = 3.500000000000000e-01f;
const real K_mNa = 8.750000000000000e+01f;
const real K_mCa = 1.380000000000000e+00f;
const real K_sat = 1.000000000000000e-01f;
const real i_CaP_max = 2.750000000000000e-01f;
const real K_rel = 3.000000000000000e+01f;
const real tau_tr = 1.800000000000000e+02f;
const real I_up_max = 5.000000000000000e-03f;
const real K_up = 9.200000000000000e-04f;
const real Ca_up_max = 1.500000000000000e+01f;
const real CMDN_max = 5.000000000000000e-02f;
const real Km_CMDN = 2.380000000000000e-03f;
const real TRPN_max = 7.000000000000001e-02f;
const real Km_TRPN = 5.000000000000000e-04f;
const real CSQN_max = 1.000000000000000e+01f;
const real Km_CSQN = 8.000000000000000e-01f;
const real V_cell = 2.010000000000000e+04f;

real calc_E_Na = (((R*T)/F)*log((Na_o/Na_i_old_)));	//3
real calc_alpha_m = 0.0f;
IFNUMBER_1(calc_alpha_m);	//4
real calc_beta_m = (8.000000000000000e-02*exp(((-V_old_)/1.100000000000000e+01)));	//5
real calc_alpha_h = 0.0f;
IFNUMBER_2(calc_alpha_h);	//9
real calc_beta_h = 0.0f;
IFNUMBER_3(calc_beta_h);	//10
real calc_alpha_j = 0.0f;
IFNUMBER_4(calc_alpha_j);	//14
real calc_beta_j = 0.0f;
IFNUMBER_5(calc_beta_j);	//15
real calc_E_K = (((R*T)/F)*log((K_o/K_i_old_)));	//19
real calc_alpha_oa = (6.500000000000000e-01*pow((exp(((V_old_-(-1.000000000000000e+01))/(-8.500000000000000e+00)))+exp((((V_old_-(-1.000000000000000e+01))-4.000000000000000e+01)/(-5.900000000000000e+01)))),(-1.000000000000000e+00)));	//22
real calc_beta_oa = (6.500000000000000e-01*pow((2.500000000000000e+00+exp((((V_old_-(-1.000000000000000e+01))+7.200000000000000e+01)/1.700000000000000e+01))),(-1.000000000000000e+00)));	//23
real calc_oa_infinity = pow((1.000000000000000e+00+exp((((V_old_-(-1.000000000000000e+01))+1.047000000000000e+01)/(-1.754000000000000e+01)))),(-1.000000000000000e+00));	//25
real calc_alpha_oi = pow((1.853000000000000e+01+(1.000000000000000e+00*exp((((V_old_-(-1.000000000000000e+01))+1.037000000000000e+02)/1.095000000000000e+01)))),(-1.000000000000000e+00));	//27
real calc_beta_oi = pow((3.556000000000000e+01+(1.000000000000000e+00*exp((((V_old_-(-1.000000000000000e+01))-8.740000000000000e+00)/(-7.440000000000000e+00))))),(-1.000000000000000e+00));	//28
real calc_oi_infinity = pow((1.000000000000000e+00+exp((((V_old_-(-1.000000000000000e+01))+3.310000000000000e+01)/5.300000000000000e+00))),(-1.000000000000000e+00));	//30
real calc_g_Kur = (5.000000000000000e-03+(5.000000000000000e-02/(1.000000000000000e+00+exp(((V_old_-1.500000000000000e+01)/(-1.300000000000000e+01))))));	//32
real calc_alpha_ua = (6.500000000000000e-01*pow((exp(((V_old_-(-1.000000000000000e+01))/(-8.500000000000000e+00)))+exp((((V_old_-(-1.000000000000000e+01))-4.000000000000000e+01)/(-5.900000000000000e+01)))),(-1.000000000000000e+00)));	//34
real calc_beta_ua = (6.500000000000000e-01*pow((2.500000000000000e+00+exp((((V_old_-(-1.000000000000000e+01))+7.200000000000000e+01)/1.700000000000000e+01))),(-1.000000000000000e+00)));	//35
real calc_ua_infinity = pow((1.000000000000000e+00+exp((((V_old_-(-1.000000000000000e+01))+2.030000000000000e+01)/(-9.600000000000000e+00)))),(-1.000000000000000e+00));	//37
real calc_alpha_ui = pow((2.100000000000000e+01+(1.000000000000000e+00*exp((((V_old_-(-1.000000000000000e+01))-1.950000000000000e+02)/(-2.800000000000000e+01))))),(-1.000000000000000e+00));	//39
real calc_beta_ui = (1.000000000000000e+00/exp((((V_old_-(-1.000000000000000e+01))-1.680000000000000e+02)/(-1.600000000000000e+01))));	//40
real calc_ui_infinity = pow((1.000000000000000e+00+exp((((V_old_-(-1.000000000000000e+01))-1.094500000000000e+02)/2.748000000000000e+01))),(-1.000000000000000e+00));	//42
real calc_alpha_xr = 0.0f;
IFNUMBER_6(calc_alpha_xr);	//45
real calc_beta_xr = 0.0f;
IFNUMBER_7(calc_beta_xr);	//46
real calc_xr_infinity = pow((1.000000000000000e+00+exp(((V_old_+1.410000000000000e+01)/(-6.500000000000000e+00)))),(-1.000000000000000e+00));	//48
real calc_alpha_xs = 0.0f;
IFNUMBER_8(calc_alpha_xs);	//51
real calc_beta_xs = 0.0f;
IFNUMBER_9(calc_beta_xs);	//52
real calc_xs_infinity = pow((1.000000000000000e+00+exp(((V_old_-1.990000000000000e+01)/(-1.270000000000000e+01)))),(-5.000000000000000e-01));	//54
real calc_i_Ca_L = (Cm*g_Ca_L*d_old_*f_old_*f_Ca_old_*(V_old_-6.500000000000000e+01));	//56
real calc_d_infinity = pow((1.000000000000000e+00+exp(((V_old_+1.000000000000000e+01)/(-8.000000000000000e+00)))),(-1.000000000000000e+00));	//57
real calc_tau_d = 0.0f;
IFNUMBER_10(calc_tau_d);	//58
real calc_f_infinity = (exp(((-(V_old_+2.800000000000000e+01))/6.900000000000000e+00))/(1.000000000000000e+00+exp(((-(V_old_+2.800000000000000e+01))/6.900000000000000e+00))));	//60
real calc_tau_f = (9.000000000000000e+00*pow(((1.970000000000000e-02*exp(((-pow(3.370000000000000e-02,2.000000000000000e+00))*pow((V_old_+1.000000000000000e+01),2.000000000000000e+00))))+2.000000000000000e-02),(-1.000000000000000e+00)));	//61
real calc_f_Ca_infinity = pow((1.000000000000000e+00+(Ca_i_old_/3.500000000000000e-04)),(-1.000000000000000e+00));	//63
real calc_tau_f_Ca = 2.000000000000000e+00;	//64
real calc_sigma = ((1.000000000000000e+00/7.000000000000000e+00)*(exp((Na_o/6.730000000000000e+01))-1.000000000000000e+00));	//66
real calc_E_Ca = (((R*T)/(2.000000000000000e+00*F))*log((Ca_o/Ca_i_old_)));	//69
real calc_i_NaCa = ((Cm*I_NaCa_max*((exp(((gamma*F*V_old_)/(R*T)))*pow(Na_i_old_,3.000000000000000e+00)*Ca_o)-(exp((((gamma-1.000000000000000e+00)*F*V_old_)/(R*T)))*pow(Na_o,3.000000000000000e+00)*Ca_i_old_)))/((pow(K_mNa,3.000000000000000e+00)+pow(Na_o,3.000000000000000e+00))*(K_mCa+Ca_o)*(1.000000000000000e+00+(K_sat*exp((((gamma-1.000000000000000e+00)*V_old_*F)/(R*T)))))));	//73
real calc_i_CaP = ((Cm*i_CaP_max*Ca_i_old_)/(5.000000000000000e-04+Ca_i_old_));	//74
real calc_i_rel = (K_rel*pow(u_old_,2.000000000000000e+00)*v_old_*w_old_*(Ca_rel_old_-Ca_i_old_));	//76
real calc_tau_u = 8.000000000000000e+00;	//77
real calc_tau_w = 0.0f;
IFNUMBER_11(calc_tau_w);	//83
real calc_w_infinity = (1.000000000000000e+00-pow((1.000000000000000e+00+exp(((-(V_old_-4.000000000000000e+01))/1.700000000000000e+01))),(-1.000000000000000e+00)));	//84
real calc_i_tr = ((Ca_up_old_-Ca_rel_old_)/tau_tr);	//86
real calc_i_up = (I_up_max/(1.000000000000000e+00+(K_up/Ca_i_old_)));	//87
real calc_i_up_leak = ((I_up_max*Ca_up_old_)/Ca_up_max);	//88
//real calc_Ca_CMDN = ((CMDN_max*Ca_i_old_)/(Ca_i_old_+Km_CMDN));	//89
//real calc_Ca_TRPN = ((TRPN_max*Ca_i_old_)/(Ca_i_old_+Km_TRPN));	//90
//real calc_Ca_CSQN = ((CSQN_max*Ca_rel_old_)/(Ca_rel_old_+Km_CSQN));	//91
real calc_V_i = (V_cell*6.800000000000000e-01);	//92
real calc_V_rel = (4.800000000000000e-03*V_cell);	//93
real calc_V_up = (5.520000000000000e-02*V_cell);	//94
real calc_B2 = (1.000000000000000e+00+((TRPN_max*Km_TRPN)/pow((Ca_i_old_+Km_TRPN),2.000000000000000e+00))+((CMDN_max*Km_CMDN)/pow((Ca_i_old_+Km_CMDN),2.000000000000000e+00)));	//99
real calc_m_inf = (calc_alpha_m/(calc_alpha_m+calc_beta_m));	//6
real calc_tau_m = (1.000000000000000e+00/(calc_alpha_m+calc_beta_m));	//7
real calc_h_inf = (calc_alpha_h/(calc_alpha_h+calc_beta_h));	//11
real calc_tau_h = (1.000000000000000e+00/(calc_alpha_h+calc_beta_h));	//12
real calc_j_inf = (calc_alpha_j/(calc_alpha_j+calc_beta_j));	//16
real calc_tau_j = (1.000000000000000e+00/(calc_alpha_j+calc_beta_j));	//17
real calc_i_K1 = ((Cm*g_K1*(V_old_-calc_E_K))/(1.000000000000000e+00+exp((7.000000000000001e-02*(V_old_+8.000000000000000e+01)))));	//20
real calc_i_to = (Cm*g_to*pow(oa_old_,3.000000000000000e+00)*oi_old_*(V_old_-calc_E_K));	//21
real calc_tau_oa = (pow((calc_alpha_oa+calc_beta_oa),(-1.000000000000000e+00))/K_Q10);	//24
real calc_tau_oi = (pow((calc_alpha_oi+calc_beta_oi),(-1.000000000000000e+00))/K_Q10);	//29
real calc_i_Kur = (Cm*calc_g_Kur*pow(ua_old_,3.000000000000000e+00)*ui_old_*(V_old_-calc_E_K));	//33
real calc_tau_ua = (pow((calc_alpha_ua+calc_beta_ua),(-1.000000000000000e+00))/K_Q10);	//36
real calc_tau_ui = (pow((calc_alpha_ui+calc_beta_ui),(-1.000000000000000e+00))/K_Q10);	//41
real calc_i_Kr = ((Cm*g_Kr*xr_old_*(V_old_-calc_E_K))/(1.000000000000000e+00+exp(((V_old_+1.500000000000000e+01)/2.240000000000000e+01))));	//44
real calc_tau_xr = pow((calc_alpha_xr+calc_beta_xr),(-1.000000000000000e+00));	//47
real calc_i_Ks = (Cm*g_Ks*pow(xs_old_,2.000000000000000e+00)*(V_old_-calc_E_K));	//50
real calc_tau_xs = (5.000000000000000e-01*pow((calc_alpha_xs+calc_beta_xs),(-1.000000000000000e+00)));	//53
real calc_f_NaK = pow((1.000000000000000e+00+(1.245000000000000e-01*exp((((-1.000000000000000e-01)*F*V_old_)/(R*T))))+(3.650000000000000e-02*calc_sigma*exp((((-F)*V_old_)/(R*T))))),(-1.000000000000000e+00));	//67
real calc_i_B_Na = (Cm*g_B_Na*(V_old_-calc_E_Na));	//70
real calc_i_B_Ca = (Cm*g_B_Ca*(V_old_-calc_E_Ca));	//71
real calc_i_B_K = (Cm*g_B_K*(V_old_-calc_E_K));	//72
real calc_i_Na = (Cm*g_Na*pow(m_old_,3.000000000000000e+00)*h_old_*j_old_*(V_old_-calc_E_Na));	//2
real calc_i_NaK = ((((Cm*i_NaK_max*calc_f_NaK*1.000000000000000e+00)/(1.000000000000000e+00+pow((Km_Na_i/Na_i_old_),1.500000000000000e+00)))*K_o)/(K_o+Km_K_o));	//68
real calc_Fn = (1.000000000000000e+03*((1.000000000000000e-15*calc_V_rel*calc_i_rel)-((1.000000000000000e-15/(2.000000000000000e+00*F))*((5.000000000000000e-01*calc_i_Ca_L)-(2.000000000000000e-01*calc_i_NaCa)))));	//75
real calc_B1 = ((((2.000000000000000e+00*calc_i_NaCa)-(calc_i_CaP+calc_i_Ca_L+calc_i_B_Ca))/(2.000000000000000e+00*calc_V_i*F))+(((calc_V_up*(calc_i_up_leak-calc_i_up))+(calc_i_rel*calc_V_rel))/calc_V_i));	//98
real calc_u_infinity = pow((1.000000000000000e+00+exp(((-(calc_Fn-3.417500000000000e-13))/1.367000000000000e-15))),(-1.000000000000000e+00));	//78
real calc_tau_v = (1.910000000000000e+00+(2.090000000000000e+00*pow((1.000000000000000e+00+exp(((-(calc_Fn-3.417500000000000e-13))/1.367000000000000e-15))),(-1.000000000000000e+00))));	//80
real calc_v_infinity = (1.000000000000000e+00-pow((1.000000000000000e+00+exp(((-(calc_Fn-6.834999999999999e-14))/1.367000000000000e-15))),(-1.000000000000000e+00)));	//81

// Full Euler
//rDY[0] = ((-(calc_i_Na+calc_i_K1+calc_i_to+calc_i_Kur+calc_i_Kr+calc_i_Ks+calc_i_B_Na+calc_i_B_Ca+calc_i_NaK+calc_i_CaP+calc_i_NaCa+calc_i_Ca_L+stim_current))/Cm);
//rDY[1] = ((calc_m_inf-m_old_)/calc_tau_m);
//rDY[2] = ((calc_h_inf-h_old_)/calc_tau_h);
//rDY[3] = ((calc_j_inf-j_old_)/calc_tau_j);
//rDY[4] = ((calc_oa_infinity-oa_old_)/calc_tau_oa);
//rDY[5] = ((calc_oi_infinity-oi_old_)/calc_tau_oi);
//rDY[6] = ((calc_ua_infinity-ua_old_)/calc_tau_ua);
//rDY[7] = ((calc_ui_infinity-ui_old_)/calc_tau_ui);
//rDY[8] = ((calc_xr_infinity-xr_old_)/calc_tau_xr);
//rDY[9] = ((calc_xs_infinity-xs_old_)/calc_tau_xs);
//rDY[10] = ((calc_d_infinity-d_old_)/calc_tau_d);
//rDY[11] = ((calc_f_infinity-f_old_)/calc_tau_f);
//rDY[12] = ((calc_f_Ca_infinity-f_Ca_old_)/calc_tau_f_Ca);
//rDY[13] = ((calc_u_infinity-u_old_)/calc_tau_u);
//rDY[14] = ((calc_v_infinity-v_old_)/calc_tau_v);
//rDY[15] = ((calc_w_infinity-w_old_)/calc_tau_w);
//rDY[16] = ((((-3.000000000000000e+00)*calc_i_NaK)-((3.000000000000000e+00*calc_i_NaCa)+calc_i_B_Na+calc_i_Na))/(calc_V_i*F));
//rDY[17] = (((2.000000000000000e+00*calc_i_NaK)-(calc_i_K1+calc_i_to+calc_i_Kur+calc_i_Kr+calc_i_Ks+calc_i_B_K))/(calc_V_i*F));
//rDY[18] = (calc_B1/calc_B2);
//rDY[19] = (calc_i_up-(calc_i_up_leak+((calc_i_tr*calc_V_rel)/calc_V_up)));
//rDY[20] = ((calc_i_tr-calc_i_rel)*pow((1.000000000000000e+00+((CSQN_max*Km_CSQN)/pow((Ca_rel_old_+Km_CSQN),2.000000000000000e+00))),(-1.000000000000000e+00)));

// Euler + Rush-Larsen
rDY[0] = ((-(calc_i_Na+calc_i_K1+calc_i_to+calc_i_Kur+calc_i_Kr+calc_i_Ks+calc_i_B_Na+calc_i_B_Ca+calc_i_NaK+calc_i_CaP+calc_i_NaCa+calc_i_Ca_L+stim_current))/Cm);

rDY[1] = calc_m_inf + (m_old_-calc_m_inf)*exp(-dt/calc_tau_m);
rDY[2] = calc_h_inf + (h_old_-calc_h_inf)*exp(-dt/calc_tau_h);
rDY[3] = calc_j_inf + (j_old_-calc_j_inf)*exp(-dt/calc_tau_j);
rDY[4] = calc_oa_infinity + (oa_old_-calc_oa_infinity)*exp(-dt/calc_tau_oa);
rDY[5] = calc_oi_infinity + (oi_old_-calc_oi_infinity)*exp(-dt/calc_tau_oi);
rDY[6] = calc_ua_infinity + (ua_old_-calc_ua_infinity)*exp(-dt/calc_tau_ua);
rDY[7] = calc_ui_infinity + (ui_old_-calc_ui_infinity)*exp(-dt/calc_tau_ui);
rDY[8] = calc_xr_infinity + (xr_old_-calc_xr_infinity)*exp(-dt/calc_tau_xr);
rDY[9] = calc_xs_infinity + (xs_old_-calc_xs_infinity)*exp(-dt/calc_tau_xs);
rDY[10] = calc_d_infinity + (d_old_-calc_d_infinity)*exp(-dt/calc_tau_d);
rDY[11] = calc_f_infinity + (f_old_-calc_f_infinity)*exp(-dt/calc_tau_f);
rDY[12] = calc_f_Ca_infinity + (f_Ca_old_-calc_f_Ca_infinity)*exp(-dt/calc_tau_f_Ca);
rDY[13] = calc_u_infinity + (u_old_-calc_u_infinity)*exp(-dt/calc_tau_u);
rDY[14] = calc_v_infinity + (v_old_-calc_v_infinity)*exp(-dt/calc_tau_v);
rDY[15] = calc_w_infinity + (w_old_-calc_w_infinity)*exp(-dt/calc_tau_w);

rDY[16] = ((((-3.000000000000000e+00)*calc_i_NaK)-((3.000000000000000e+00*calc_i_NaCa)+calc_i_B_Na+calc_i_Na))/(calc_V_i*F));
rDY[17] = (((2.000000000000000e+00*calc_i_NaK)-(calc_i_K1+calc_i_to+calc_i_Kur+calc_i_Kr+calc_i_Ks+calc_i_B_K))/(calc_V_i*F));
rDY[18] = (calc_B1/calc_B2);
rDY[19] = (calc_i_up-(calc_i_up_leak+((calc_i_tr*calc_V_rel)/calc_V_up)));
rDY[20] = ((calc_i_tr-calc_i_rel)*pow((1.000000000000000e+00+((CSQN_max*Km_CSQN)/pow((Ca_rel_old_+Km_CSQN),2.000000000000000e+00))),(-1.000000000000000e+00)));
