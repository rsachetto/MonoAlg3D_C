//Drug
real ACh         = 0.0;
real IK1_Block   = 0.0;
real IKr_Block   = 0.0;
real IKs_Block   = 0.0;
real IKur_Block  = 0.0;
real INa_Block   = 0.0;
real INaK_Block  = 0.0;
real INCX_Block  = 0.0;
real Ito_Block   = 0.0;
real ICaL_Block  = 0.0;
real IK2P_Block  = 0.0;
real IKACH_Block = 0.0;

// POM-specific
real GNa  = extra_parameters[ 8] * (1 - INa_Block);          
real Gto  = extra_parameters[ 2] * (1 - Ito_Block); 
real GCaL = extra_parameters[ 5] * (1 - ICaL_Block);    
real GKur = extra_parameters[ 0] * (1 - IKur_Block);  
real GKr  = extra_parameters[ 1] * (1 - IKr_Block);           
real GKs  = extra_parameters[ 4] * (1 - IKs_Block);           
real GK1  = extra_parameters[ 3] * (1 - IK1_Block);
real GNCX = extra_parameters[ 7] * (1 - INCX_Block);          
real GCaP = 1.0;
real GNab = 1.0;
real GCab = 1.0;
real GNaK = extra_parameters[ 6] * (1 - INaK_Block);
real GK2P = extra_parameters[ 9] * (1 - IK2P_Block);
real GKACH= extra_parameters[10] * (1 - IKACH_Block);          
real Grel = 1.0; 
real Gup  = 1.0; 
real Gleak= 1.0;

// Region-specific
// Right Atrium
if        (mapping == 1){
    GKr                     = GKr  * (1/1.6);
// Crista Terminalis
} else if (mapping == 2) {
    Gto                     = Gto  * 1.35;
    GCaL                    = GCaL * 1.60;
    GKr                     = GKr  * 0.90;
// Pectinate Muscles
} else if (mapping == 3) {
    Gto                     = Gto  * 1.05;
    GCaL                    = GCaL * 0.95;
    GKr                     = GKr  * 0.90;
// Tricuspid annulus
} else if (mapping == 5) {
    Gto                     = Gto  * 1.05;
    GCaL                    = GCaL * 0.65;
    GKr                     = GKr  * 3.00;
// Left atrial appendage
} else if (mapping == 8) {
    Gto                     = Gto  * 0.65;
    GCaL                    = GCaL * 1.05;
    GKr                     = GKr  * 2.75; 
// Mitral annulus
} else if (mapping == 9) {
    Gto                     = Gto  * 1.05;
    GCaL                    = GCaL * 0.65;
    GKr                     = GKr  * 3.00;    


// Sinoatrial node + GP
} else if (mapping == 100){
    ACh                     = 0.1;
// Right Atrium + GP
} else if (mapping == 101){
    GKr                     = GKr  * (1/1.6);
    ACh                     = 0.1;
// Crista Terminalis + GP
} else if (mapping == 102) {
    Gto                     = Gto  * 1.35;
    GCaL                    = GCaL * 1.60;
    GKr                     = GKr  * 0.90;
    ACh                     = 0.1;
// Pectinate Muscles + GP
} else if (mapping == 103) {
    Gto                     = Gto  * 1.05;
    GCaL                    = GCaL * 0.95;
    GKr                     = GKr  * 0.90;
    ACh                     = 0.1;
// Right atrial appendage + GP
} else if (mapping == 104) {
    ACh                     = 0.1;
// Tricuspid annulus + GP
} else if (mapping == 105) {
    Gto                     = Gto  * 1.05;
    GCaL                    = GCaL * 0.65;
    GKr                     = GKr  * 3.00;
    ACh                     = 0.1;
// Right atrial appendage + GP
} else if (mapping == 106) {
    ACh                     = 0.1;
// Right atrial appendage + GP
} else if (mapping == 107) {
    ACh                     = 0.1;
// Left atrial appendage + GP
} else if (mapping == 108) {
    Gto                     = Gto  * 0.65;
    GCaL                    = GCaL * 1.05;
    GKr                     = GKr  * 2.75;
    ACh                     = 0.1; 
// Mitral annulus
} else if (mapping == 109) {
    Gto                     = Gto  * 1.05;
    GCaL                    = GCaL * 0.65;
    GKr                     = GKr  * 3.00;
    ACh                     = 0.1;  
// Mitral annulus
} else if (mapping == 110) {
    Gto                     = Gto  * 1.05;
    GCaL                    = GCaL * 0.65;
    GKr                     = GKr  * 3.00;
    ACh                     = 0.1;     
}

// Parameters 
const real F                 	= 96.4867;
const real T                 	= 310;                    // K (37 C)
const real R                  	= 8.31441;                // J mol^-1 K^-1
const real FoRT               	= F/(R*T);

// Capacitance and cell structure
const real Cm                	= 100;                    // pF
const real Vcell            	= 20100.0;                // um^3
const real Vcyto              	= 0.68   * Vcell;         // um^3
const real VjSR               	= 0.0048 * Vcell;         // um^3
const real VnSR                 = 0.0552 * Vcell;         // um^3

// Current parameters
const real gNa               	= 7.8;                    // s/mF
const real gto              	= 0.1652;                 // s/mF
const real gCaL                 = 0.1238;                 // s/mF 
const real gKr                  = 0.0294*sqrt(Ko/5.4); 	  // s/mF
const real gKs               	= 0.129;                  // s/mF
const real gK1              	= 0.09*pow(Ko/5.4,0.4);   // s/mF
const real gK2P                 = 0.0035;                 // s/mF
    
const real INCX_bar             = 1600;                   // pA/pF
const real INCX_kNao        	= 87.5;                   // mM
const real INCX_kCao         	= 1.38;                   // mM
const real INCX_k            	= 0.1;
const real INCX_gamma         	= 0.35;

const real INaK_bar         	= 0.59933874;             // pA/pF
const real INaK_kK              = 1.5;                    // mM
const real INaK_kNa           	= 10.0;                   // mM

const real ICaP_bar             = 0.275;                  // pA/pF
const real ICaP_kCa           	= 0.5e-3;                 // mM
const real gNab              	= 0.000674;               // s/mF
const real gCab             	= 0.001131;               // s/mF

const real ICaL_ci_tau        	= 2.0;                    // ms

// Ca2+ handling
const real cmdnbar              = 0.050;                  // mM
const real trpnbar           	= 0.070;                  // mM
const real csqnbar           	= 10.0;		              // mM
const real cmdn_k            	= 0.00238;                // mM
const real trpn_k             	= 0.0005;                 // mM

const real J_rel_max            = 30.0;                   // mM/ms
const real J_SERCA_max       	= 0.005;                  // mM/ms
const real J_SERCA_kCa       	= 0.00092;                // mM
const real J_leak_max        	= 0.005;                  // mM/ms
const real J_leak_kCaSR       	= 15.0;                   // mM
const real J_jsr_nsr_tau     	= 180.0;                  // ms

// Reversal Potentials
real ENa                        = 	((R * T)/F)*log(Nao/Nai);
real EK                         = 	((R * T)/F)*log(Ko /Ki );
real ECa                        = 0.5 *	((R * T)/F)*log(Cao/Cai);

// INa Luo Rudy
// Set Activation gate alpha and beta
real INa_va_al                  = (fabs(Vm+47.13)<1e-10 ? 3.2 : 0.32*(Vm+47.13)/(1-exp(-0.1*(Vm+47.13)))) ;
real INa_va_bet                 = 0.08*exp(-Vm/11.0);
// Set inactivation gates alphas and betas
real INa_vi_1_al                = (Vm < -40.0 ? 0.135*exp((80+Vm)/-6.8) : 0);
real INa_vi_1_bet               = (Vm < -40.0 ? 3.56*exp(0.079*Vm)+310000*exp(0.35*Vm) : 1.0/(0.13*(1+exp((Vm+10.66)/-11.1))));
real INa_vi_2_al                = (Vm < -40.0 ? (-127140*exp(0.2444*Vm)-0.00003474*exp(-0.04391*Vm))*((Vm+37.78)/(1+exp(0.311*(Vm+79.23)))) : 0);
real INa_vi_2_bet               = (Vm < -40.0 ? (0.1212*exp(-0.01052*Vm))/(1+exp(-0.1378*(Vm+40.14))) : (0.3*exp(-0.0000002535*Vm))/(1+exp(-0.1*(Vm+32))));
// Set tau and SS from alpha and beta
real INa_va_tau      	        = 1.0/(INa_va_al + INa_va_bet);     // 1/(a+b)
real INa_vi_1_tau      	        = 1.0/(INa_vi_1_al + INa_vi_1_bet);
real INa_vi_2_tau               = 1.0/(INa_vi_2_al + INa_vi_2_bet);
real INa_va_ss         	        = INa_va_al   * INa_va_tau;           // a*tau
real INa_vi_1_ss                = INa_vi_1_al * INa_vi_1_tau;
real INa_vi_2_ss                = INa_vi_2_al * INa_vi_2_tau;
// Current
real INa                        = GNa * gNa * pow(INa_va, 3) * INa_vi_1 * INa_vi_2 * (Vm - ENa);

// Ito 
// Voltage activation	
real Ito_va_ss                  = 1/(1+exp(-(Vm+20.47)/17.54));
real Ito_va_al         	        = 0.65/(exp(-(Vm+10)/8.5)+exp(-(Vm-30)/59));
real Ito_va_bet                 = 0.65/(2.5+exp((Vm+82)/17));
real Ito_va_tau                 = 1.0/(3*(Ito_va_al+Ito_va_bet));
// Voltage inactivation
real Ito_vi_ss                  = 1/(1+exp((Vm+43.1)/5.3));
real Ito_vi_al                  = 1.0/(18.53+exp((Vm+113.7)/10.95));
real Ito_vi_bet                 = 1.0/(35.56+exp(-(Vm+1.26)/7.44));
real Ito_vi_tau                 = 1.0/(3*(Ito_vi_al +  Ito_vi_bet)); 
// Current
real Ito                        = Gto * gto * pow(Ito_va, 3) * Ito_vi * (Vm - EK);

// ICaL
real ICaL_ci_ss                 = 1/(1+Cai/0.00035);
// Voltage activation 
real ICaL_va_ss                 = 1/(1+exp((Vm+10)/-8));
real ICaL_va_al                 = 1/(1+exp((Vm+10)/-6.24));
real ICaL_va_be                 =   (1-exp((Vm+10)/-6.24))/(0.035*(Vm+10));
real ICaL_va_tau 				= (fabs(Vm+10)<1.0e-10 ?  4.579*ICaL_va_al : ICaL_va_al*ICaL_va_be);
// Voltage inactivation 
real ICaL_vi_ss                 = exp(-(Vm+28.0)/6.9)/(1.0+exp(-(Vm+28.0)/6.9));
real ICaL_vi_tau                = 9.0/(0.0197*exp(-pow(0.0337,2)*pow((Vm+10),2))+0.02);
// Current
real ICaL                       = GCaL * gCaL * ICaL_va * ICaL_vi * ICaL_ci * (Vm - 65);

// IKur
// Voltage activation   
real IKur_va_ss                 = 1/(1+exp(-(Vm+30.3)/9.6));
real IKur_va_al                 = 0.65/(exp(-(Vm+10)/8.5)+exp(-(Vm-30)/59.0));
real IKur_va_bet                = 0.65/(2.5+exp((Vm+82)/17.0));
real IKur_va_tau                = 1.0/(3*(IKur_va_al + IKur_va_bet));
// Voltage inactivation
real IKur_viS_ss                = 1/(1+exp((Vm+5)/5));
real IKur_viS_tau               = 5800 * (1+exp(-(Vm+80)/11));
real IKur_viF_ss                = 1/(1+exp((Vm-35)/20));
real IKur_viF_tau               = 800 * (2 - Vm/40);
// Current
real IKur_dynamic_g             = 0.005+0.05/(1+exp(-(Vm-15)/13));
real IKur                       = GKur * IKur_dynamic_g * pow(IKur_va, 3) * IKur_viS * IKur_viF * (Vm - EK);

// IK1
real IK1_va_ti                  = (1.0 + exp(0.07*(Vm+80)));
real IK1                        = GK1 * gK1 * (Vm - EK)/IK1_va_ti;

// IKr
// Voltage activation  
real IKr_va_ss                  = 1/(1+exp(-(Vm+14.1)/6.5));
real IKr_va_al                  = 0.0003      * (Vm+14.1)/(1-exp(-(Vm+14.1)/5));
real IKr_va_be                  = 0.000073898 * (Vm-3.3328)/(exp((Vm-3.3328)/5.1237)-1);
real IKr_va_tau                 = 1/(IKr_va_al+IKr_va_be);
// Time-independant gate
real IKr_vi_ti                  = 1/(1+exp((Vm+15)/22.4));
real IKr                        = GKr * gKr * IKr_va *  IKr_vi_ti * (Vm - EK);

// IKs
// Voltage activation   
real IKs_va_ss                  = 1.0/pow((1+exp(-(Vm-19.9)/12.7)),0.5);
real IKs_va_al                  = 0.00004  * (Vm-19.9)/(1-exp(-(Vm-19.9)/17));
real IKs_va_be                  = 0.000035 * (Vm-19.9)/(exp((Vm-19.9)/9)-1);
real IKs_va_tau                 = 0.5/(IKs_va_al+IKs_va_be);
// Current
real IKs                        = GKs * gKs * IKs_va * IKs_va * (Vm - EK);

// INCX
real INCX                       = GNCX * INCX_bar *( exp(INCX_gamma* FoRT*Vm) * Nai*Nai*Nai * Cao - exp((INCX_gamma-1)* FoRT*Vm) * Nao*Nao*Nao*Cai)/((pow(INCX_kNao,3) + Nao*Nao*Nao)*(INCX_kCao + Cao) * (1 + INCX_k * exp((INCX_gamma-1)* FoRT*Vm)));

// INaK
real sigma                      = (exp(Nao/67.3)-1.0)/7.0;
real FNaK                       = pow(1.0+0.1245*exp(-0.1*Vm*FoRT)+0.0365*sigma*exp(-Vm*FoRT ), -1.0);
real INaK                       = GNaK * INaK_bar * FNaK * (1.0/(1.0+pow((INaK_kNa /Nai),1.5))) * (Ko/(Ko + INaK_kK));

// ICap
real ICaP                       = GCaP * (ICaP_bar * Cai) / (Cai + ICaP_kCa);

// INab
real INab                       = GNab * gNab * (Vm - ENa);

// ICab
real ICab                       = GCab * gCab * (Vm - ECa);

// IK2P
// Voltage activation
real IK2P_va_ss                 = 0.2 + 0.8/(1+exp(-(Vm+5)/14));
real IK2P_va_tau                = 2.0 +  40/(1+exp(pow(Vm+25,2)/80));
// Current
real IK2P                       = GK2P * gK2P * IK2P_va * (Vm - EK);

// IKACH
real Gkach                      = GKACH;
real r_kach                     = 0.055+5.0/(1+exp((Vm+85.0)/5.0));
real d_kach                     = 10.0/(1+9.14/pow(ACh*0.1,0.478));
// current
real IKACh                      = Gkach*d_kach*r_kach*(Vm-EK);


// Stimulation
real Istim                      = stim_current;

// Itot
real Itot                       = INa + Ito + IK1 + ICaL + IKur + INCX + INaK + ICaP + INab +ICab + IKr + IKs + IK2P + IKACh;


// Ca Homeostasis
// Jrel
real J_rel                      = Grel * J_rel_max * RyRo*RyRo * RyRr * RyRi *(CajSR - Cai);
real Cai_Fn                     = 1.0e-12*VjSR*J_rel - 1.0e-12/(2.0*F)*(0.5*ICaL-0.2*INCX)*Cm;

// RyR
real RyRo_ss                    =     1.0/(1.0+exp(-(Cai_Fn-3.4175e-13)/13.67e-16));
real RyRr_ss                    = 1.0-1.0/(1.0+exp(-(Cai_Fn-6.835e-14 )/13.67e-16));
real RyRi_ss                    = 1.0-1.0/(1.0+exp(-(Vm-40.0)/17.0));

real RyRo_tau                   = 8.0;
real RyRr_tau                   = 1.91+2.09/(1.0+exp(-(Cai_Fn-3.4175e-13)/13.67e-16));
real RyRi_tau                   = 6.0*(1.0-exp(-(Vm -7.9)/5.0))/((1.0+0.3*exp(-(Vm -7.9)/5.0))*(Vm -7.9));

// Jup
real J_SERCA                    = Gup   * J_SERCA_max * Cai/(Cai + J_SERCA_kCa);
real J_leak                     = Gleak * J_leak_max  * CanSR / J_leak_kCaSR;

// SR compart transfer
real J_jsr_nsr                  = (CanSR - CajSR)/J_jsr_nsr_tau;

// UPDATE
rDY[0]  = -(Itot+Istim);
// INa
rDY[1]  = (INa_va_ss   - INa_va)   / INa_va_tau;
rDY[2]  = (INa_vi_1_ss - INa_vi_1) / INa_vi_1_tau;
rDY[3]  = (INa_vi_2_ss - INa_vi_2) / INa_vi_2_tau;
// Ito
rDY[4]  = (Ito_va_ss   - Ito_va)   / Ito_va_tau;
rDY[5]  = (Ito_vi_ss   - Ito_vi)   / Ito_vi_tau;
// ICaL
rDY[6]  = (ICaL_va_ss  - ICaL_va)  / ICaL_va_tau;
rDY[7]  = (ICaL_vi_ss  - ICaL_vi)  / ICaL_vi_tau;
rDY[8]  = (ICaL_ci_ss  - ICaL_ci)  / ICaL_ci_tau;
// IKur
rDY[9]  = (IKur_va_ss  - IKur_va)  / IKur_va_tau;
rDY[10] = (IKur_viS_ss - IKur_viS) / IKur_viS_tau;
rDY[11] = (IKur_viF_ss - IKur_viF) / IKur_viF_tau;
// IKr
rDY[12] = (IKr_va_ss   - IKr_va)   / IKr_va_tau;
// IKs
rDY[13] = (IKs_va_ss   - IKs_va)   / IKs_va_tau;
// IK2P
rDY[14] = (IK2P_va_ss  - IK2P_va)  / IK2P_va_tau;
// RyR
rDY[17] = (RyRo_ss     - RyRo)     / RyRo_tau;
rDY[18] = (RyRr_ss     - RyRr)     / RyRr_tau;
rDY[19] = (RyRi_ss     - RyRi)     / RyRi_tau;
// Ca-Handling
rDY[15] = (J_jsr_nsr - J_rel)/((1+csqnbar*0.8/pow((CajSR+0.8),2)));
rDY[16] = (J_SERCA - J_jsr_nsr * VjSR/VnSR - J_leak);
// Sodium Concentratiuon
rDY[20] = Cm*(-3*INaK - 3*INCX - INab - INa)/(F*Vcyto);
rDY[21] = 0.0;
// Potassium Concentration
rDY[22] = Cm*(2*INaK - IK1 - Ito - IKur - IKr - IKs - IK2P - IKACh)/(F*Vcyto);
rDY[23] = 0.0;
// Calcium concentration
rDY[24] = Cm*(2*INCX-(ICaP+ICaL+ICab))/(2.0*Vcyto*F)+(VnSR*(J_leak - J_SERCA)+J_rel*VjSR)/Vcyto / (1.0+trpnbar*trpn_k/pow(Cai + trpn_k, 2.0) + cmdnbar*cmdn_k/pow(Cai + cmdn_k, 2.0));
rDY[25] = 0.0;
