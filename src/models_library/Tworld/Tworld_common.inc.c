/*                                      
real membrane_v;
real ions_na_junc;
real ions_na_sl;
real ions_na_i;
real ions_k_i;
real ions_ca_junc;
real ions_ca_sl;
real ions_ca_i;
real ions_cl_i;
real ions_ca_SR;
real buffers_NaBj;
real buffers_NaBsl;
real buffers_TnClow;
real buffers_TnCHc;
real buffers_TnCHm;
real buffers_CaM;
real buffers_Myosin_ca;
real buffers_Myosin_mg;
real buffers_SRB;
real buffers_SLLj;
real buffers_SLLsl;
real buffers_SLHj;
real buffers_SLHsl;
real buffers_Csqn;
real ina_m;
real ina_h;
real ina_j;
real inal_ml;
real inal_hl;
real inal_hl_p;
real ito_xtos;
real ito_ytos;
real ito_xtof;
real ito_ytof;
real ito_xtos_p;
real ito_ytos_p;
real ito_xtof_p;
real ito_ytof_p;
real iks_xs_junc;
real iks_xs_sl;
real ikr_c0;
real ikr_c1;
real ikr_c2;
real ikr_o;
real ikr_i;
real ical_d;
real ical_ff;
real ical_fs;
real ical_fcaf;
real ical_fcas;
real ical_jca;
real ical_nca;
real ical_nca_i;
real ical_ffp;
real ical_fcafp;
real ical_pureCDI_junc;
real ical_pureCDI_sl;
real ryr_R;
real ryr_O;
real ryr_I;
real ryr_CaRI;
real ryr_R_p;
real ryr_O_p;
real ryr_I_p;
real ryr_CaRI_p;
real jrel_icaldep_act;
real jrel_icaldep_f1;
real jrel_icaldep_f2;
real contraction_XS;
real contraction_XW;
real contraction_Ca_TRPN;
real contraction_TmBlocked;
real contraction_ZETAS;
real contraction_ZETAW;
real camk_trap;
real camk_f_ICaL;
real camk_f_RyR;
real camk_f_PLB;
real casig_serca_trap;
real ina_hp;
real ina_jp;
real ina_m_PKA;
real ina_h_PKA;
real ina_j_PKA;
real ina_h_dualP;
real ina_j_dualP;
real ical_d_P;
real ical_ff_P;
real ical_fs_P;
real ical_fcaf_P;
real ical_fcas_P;
real ical_fBPf;
real ical_fcaBPf;
*/



// TODO: Current multipliers (this will become extra_data later!)
real INa_Multiplier = 1.0;
real INaL_Multiplier = 1.0;
real INab_Multiplier = 1.0;
real INaK_Multiplier = 1.0;
real IKr_Multiplier = 1.0;
real IKs_Multiplier = 1.0;
real IKb_Multiplier = 1.0;
real Itos_Multiplier = 1.0;
real Itof_Multiplier = 1.0;
real IK1_Multiplier = 1.0;
real IClCa_Multiplier = 1.0;
real IClb_Multiplier = 1.0;
real ICaL_PCaMultiplier = 1.0;
real IpCa_Multiplier = 1.0;
real ICab_Multiplier = 1.0;
real INaCa_Multiplier = 1.0;
real Jrel_Multiplier = 1.0;
real Jup_Multiplier = 1.0;
real CMDN_Multiplier = 1.0;
real Whole_cell_PP1 = 0.13698;

// Land-Niederer parameters
//      input coded here:
int mode = 0;   // 0 = "intact", 1 = "skinned"
real lambda=1;
real lambda_rate=0;

int runSignallingPathway = 0;
//TODO:  PKA signalling - simulating or using direct phosphorylation values
//  Include some parameters as 'extra_data'
//  Talk to Jakub about these parameters to better understand how they work ...
real fINaP= 0;
real fICaL_PKA = 0;
real fINaK_PKA = 0;
real fIKs_PKA = 0;
real fPLB_PKA = 0;
real fTnI_PKA = 0;
real fMyBPC_PKA = 0;
if (runSignallingPathway == 1) {
    //ydot_signalling = getPKASignalling(y(93+(1:57)), constantsSig);
    //[fICaL_PKA,fIKs_PKA,fPLB_PKA,fTnI_PKA,fINa_PKA,fINaK_PKA,fRyR_PKA,fIKb_PKA]=  getEffectiveFraction(y(93+(1:57)), constantsSig); 
    //fMyBPC_PKA = fTnI_PKA;
    //// Concentration of uninhibited PP1 in the cytosolic compartment
    //pp1_PP1f_cyt_sum = constantsSig(38) - constantsSig(37) + y(93+39);
    //// Concentration of uninhibited PP1 in the cytosolic compartment
    //PP1f_cyt = 0.5 * (sqrt(pp1_PP1f_cyt_sum ^ 2.0 + 4.0 * constantsSig(38) * constantsSig(37)) - pp1_PP1f_cyt_sum);
    //Whole_cell_PP1  = constantsSig(36) / constantsSig(5) + constantsSig(35) / constantsSig(6) + PP1f_cyt / constantsSig(7);
}
    
// Default parameters
real ions_na_o = 140;
real ions_ca_o = 1.8;
real ions_k_o = 5;


real ko_NZ;
real ko_ICZ;
real fKatp_NZ;
real fKatp_ICZ;
real GNa_GCaL_NZ;
real GNa_GCaL_ICZ;

switch (NZ_severity) {
  case 0:
    ko_NZ 		= 5.0;
	fKatp_NZ 	= 0.0;
	GNa_GCaL_NZ = 1.0;
    break;
  case 1:
    ko_NZ 		= 8.0;
	fKatp_NZ 	= 0.04;
	GNa_GCaL_NZ = 0.8;
    break;
  case 2:
    ko_NZ 		= 7.0;
	fKatp_NZ 	= 0.02;
	GNa_GCaL_NZ = 0.9;
    break;
}

switch (ICZ_severity) {
  case 0:
    ko_ICZ 			= 5.0;
	fKatp_ICZ 		= 0.0;
	GNa_GCaL_ICZ 	= 1.0;
    break;
  case 1:
	ko_ICZ 			= 8.0;
	fKatp_ICZ 		= 0.04;
	GNa_GCaL_ICZ 	= 0.8;
    break;
  case 2:
	ko_ICZ 			= 7.0;
	fKatp_ICZ 		= 0.02;
	GNa_GCaL_ICZ 	= 0.9;
    break;
}


/////////////// [K+]0 /////////////////////////////////////////////////////////
	if(ischemia == 0) 	// NZ
	{
		ions_k_o = ko_NZ;
	} else
	
		if(ischemia >= 1) 	// ICZ
		{
			ions_k_o = ko_ICZ;
		} else
			// BZ, gradient
			{
				ions_k_o = ko_NZ + ((ischemia) * (ko_ICZ-ko_NZ));
			}
			


////////////////////////////////////////////////////////////////////////////////

// Localization of ICaL and NCX: the fraction in junctional subspace
real ICaL_fractionSS = 0.8; 
real INaCa_fractionSS = 0.31;

// Constants
real R = 8314;       // [J/kmol*K]
real Frdy = 96485;   // [C/mol]
real Temp = 310;     // [K]
real FoRT = (Frdy / R) / Temp;
real Cmem = 1.3810e-10;   // [F] membrane capacitance
real Qpow = (Temp-310)/10;
real F = Frdy;
real T = Temp;

// Cell geometry
real cellLength = 100;     // cell length [um]
real cellRadius = 10.25;   // cell radius [um]
real Vcell = M_PI*pow(cellRadius,2)*cellLength*1e-15;    // [L]
real Vmyo = 0.65*Vcell; 
real Vsr = 0.035*Vcell; 
real Vsl = 0.02*Vcell; 
real Vjunc = 0.0539*.01*Vcell;

real J_ca_juncsl = 1/3.06685e12;  
real J_ca_slmyo =  1/0.74556e11; 
real J_na_juncsl = 1/(1.6382e12/3*100); 
real J_na_slmyo = 1/(1.8308e10/3*100);  

// Fractional currents in compartments
real Fjunc = 0.11; 
real Fsl = 1-Fjunc;
// different localizations used for ICaL (ICaL_fractionSS), NCX (INaCa_fractionSS), and I(Ca)Cl (50:50)

// Fixed ion concentrations
real ions_cl_o = 150;  // Extracellular Cl  [mM]
real ions_mg_i = 0.5;  // Intracellular Mg  [mM]

// Nernst Potentials
real ena_junc = (1/FoRT)*log(ions_na_o/ions_na_junc);     // [mV]
real ena_sl = (1/FoRT)*log(ions_na_o/ions_na_sl);       // [mV]
real ek = (1/FoRT)*log(ions_k_o/ions_k_i);	        // [mV]
real eca_junc = (1/FoRT/2)*log(ions_ca_o/ions_ca_junc);   // [mV]
real eca_sl = (1/FoRT/2)*log(ions_ca_o/ions_ca_sl);     // [mV]
real ecl = (1/FoRT)*log(ions_cl_i/ions_cl_o);            // [mV]

//// Buffering parameters
real Bmax_Naj = 7.561;       // [mM] // Bmax_Naj = 3.7; (c-code difference?)  // Na buffering
real Bmax_Nasl = 1.65;       // [mM]
real koff_na = 1e-3;         // [1/ms]
real kon_na = 0.1e-3;        // [1/mM/ms]
real Bmax_TnClow = 70e-3;    // [mM]                      // TnC low affinity
// koff_tncl = 19.6e-3;    // [1/ms]
// kon_tncl = 32.7;        // [1/mM/ms]
real Bmax_TnChigh = 140e-3;  // [mM]                      // TnC high affinity
real koff_tnchca = 0.032e-3; // [1/ms]
real kon_tnchca = 2.37;      // [1/mM/ms]
real koff_tnchmg = 3.33e-3;  // [1/ms]
real kon_tnchmg = 3e-3;      // [1/mM/ms]
real Bmax_CaM = 24e-3 * CMDN_Multiplier;       // [mM] **? about setting to 0 in c-code**   // CaM buffering
real koff_cam = 238e-3;      // [1/ms]
real kon_cam = 34;           // [1/mM/ms]
real Bmax_myosin = 140e-3;   // [mM]                      // Myosin buffering
real koff_myoca = 0.46e-3;   // [1/ms]
real kon_myoca = 13.8;       // [1/mM/ms]
real koff_myomg = 0.057e-3;  // [1/ms]
real kon_myomg = 0.0157;     // [1/mM/ms]
real Bmax_SR = 17.85854e-3;     // [mM] (Bers text says 47e-3) 19e-3
real koff_sr = 60e-3;        // [1/ms]
real kon_sr = 100;           // [1/mM/ms]
real Bmax_SLlowsl = 33.923e-3*Vmyo/Vsl;        // [mM]    // SL buffering
real Bmax_SLlowj = 4.89983e-4*Vmyo/Vjunc;    // [mM]    //Fei *0.1!!! junction reduction factor
real koff_sll = 1300e-3;     // [1/ms]
real kon_sll = 100;          // [1/mM/ms]
real Bmax_SLhighsl = 12.15423e-3*Vmyo/Vsl;       // [mM]
real Bmax_SLhighj = 1.75755e-4*Vmyo/Vjunc;  // [mM] //Fei *0.1!!! junction reduction factor
real koff_slh = 30e-3;       // [1/ms]
real kon_slh = 100;          // [1/mM/ms]
real Bmax_Csqn = 136.55214e-3*Vmyo/Vsr;            // [mM] // Bmax_Csqn = 2.6;      // Csqn buffering
real koff_csqn = 65;         // [1/ms]
real kon_csqn = 100;         // [1/mM/ms]

// CaMK and Ca signalling
real PP1_tot = Whole_cell_PP1;       //0.13698; %Can be linked to BARS signaling later currently steady state value in the absence of ISO
// PP1_tot = pp1.PP1_cav / cell.vr_cav + pp1.PP1_eca / cell.vr_eca + pp1.PP1f_cyt / cell.vr_cyt
real CaMK0  = 2*0.05; // Equilibrium fraction of active CaMKII binding sites
real Km_CaMK_Ca = 5*0.0015; //[mmol/L] CaMKII affinity for Ca2+/CaM activation %Adjusted because of smaller cleft space

real bound = CaMK0 * (1 - camk_trap) / (1 + Km_CaMK_Ca / ions_ca_junc);
real CaMK_active = bound + camk_trap; // Fraction of active CaMKII

real alpha_camk = 0.05;
real beta_camk  = 6.8e-4;
real d_camk_trap = alpha_camk * bound * CaMK_active - beta_camk * camk_trap * (0.1 + 0.9 * PP1_tot / 0.1371); // Euler

real tau_plb = 100000; // [ms] Time constant of CaMKII PLB phosphorylation
real tau_ryr = 10000; // [ms] Time constant of CaMKII RyR phosphorylation
real tau_cal = tau_ryr; // Time constant of CaMKII ICaL phosphorylation


real K_Phos_CaMK = 0.35;  // Affinity of PLB, ICaL etc for CaMKII
real CaMK_Phos_ss_ICaL = CaMK_active / (CaMK_active + K_Phos_CaMK);
real CaMK_Phos_ss_RyR = CaMK_active / (CaMK_active + 1);
real CaMK_Phos_ss_PLB = CaMK_active / (CaMK_active + 10);

real d_camk_f_ICaL = (CaMK_Phos_ss_ICaL - camk_f_ICaL) / tau_cal;   // Rush-Larsen
real d_camk_f_RyR = (CaMK_Phos_ss_RyR - camk_f_RyR)  / tau_ryr;     // Rush-Larsen
real d_camk_f_PLB = (CaMK_Phos_ss_PLB - camk_f_PLB)  / tau_plb;     // Rush-Larsen

real alpha_serca = 0.05;
real bound_serca = CaMK0 * (1 - casig_serca_trap) / (1 + Km_CaMK_Ca / ions_ca_junc);
real casig_SERCA_act = bound_serca + casig_serca_trap; // Fraction of active CaMKII
real d_casig_serca_trap = alpha_serca * bound_serca * casig_SERCA_act - beta_camk * casig_serca_trap * (0.1 + 0.9 * PP1_tot / 0.1371); // Euler

// INa
// INa formulations - taken from Doste et al.
// The Grandi implementation updated with INa phosphorylation.
// m gate
real mss = 1 / (pow(1 + exp( -(56.86 + v) / 9.03 ),2));
real taum = 0.1292 * exp(-pow((v+45.79)/15.54,2)) + 0.06487 * exp(-pow((v-4.823)/51.12,2));
real dm = (mss - m) / taum;                     // Rush-Larsen

// h gate
real ah = (v >= -40) ? (0) : (0.057 * exp( -(v + 80) / 6.8 ));
real bh = (v >= -40) ? (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 )))) : ((2.7 * exp( 0.079 * v) + 3.1*pow(10,5) * exp(0.3485 * v)));
real tauh = 1 / (ah + bh);
real hss = 1 / (pow(1 + exp( (v + 71.55)/7.43 ),2));
real dh = (hss - h) / tauh;                     // Rush-Larsen

// j gate
real aj = (v >= -40) ? (0) : (((-2.5428 * pow(10,4)*exp(0.2444*v) - 6.948*pow(10,-6) * exp(-0.04391*v)) * (v + 37.78)) / (1 + exp( 0.311 * (v + 79.23) )));
real bj = (v >= -40) ? ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) ))) : ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) )));
real tauj = 1 / (aj + bj);
real jss = 1 / pow((1 + exp( (v + 71.55)/7.43 )),2);
real dj = (jss - j) / tauj;                     // Rush-Larsen

// gating CaMK-P
real hssp = 1 / pow((1 + exp( (v + 71.55 + 6)/7.43 )),2);
real dhp = (hssp - hp) / tauh;                  // Rush-Larsen
real taujp = 1.46 * tauj;
real djp = (jss - jp) / taujp;                  // Rush-Larsen

// gating PKA
real mss_P = 1 / (pow(1 + exp( -(56.86 + v + 5.0) / 9.03 ),2));

real dm_P = (mss_P - m_P) / taum;                     // Rush-Larsen

real hss_P = 1 / (pow(1 + exp( (v + 71.55+5.0)/7.43 ),2));  // BetaAdrenergic;
real dh_P=(hss_P-h_P)/tauh;   // BetaAdrenergic;
real jss_P=hss_P;
real dj_P = (jss_P - j_P) / tauj;               // Rush-Larsen

// Both Phosphorylated
real hssp_P = 1 / (pow(1 + exp( (v + 71.55 + 6 + 5.0)/7.43 ),2)); // BetaAdrenergic;
real dhp_P = (hssp_P - hp_P) / tauh;                 // Rush-Larsen
real jssp_P=hssp_P;
real djp_P = (jssp_P - jp_P) / taujp;                // Rush-Larsen

// Putting together the channels behavior and fraction
real GNa = 22.08788 * INa_Multiplier; 

//////////////// INa /////////////
	if(ischemia == 0) 	// NZ
	{
		GNa = GNa * GNa_GCaL_NZ;
	} else
	
		if(ischemia >= 1) 	// ICZ
		{
			GNa = GNa * GNa_GCaL_ICZ;
		} else
			// BZ, gradient
			{
				GNa = GNa * (GNa_GCaL_NZ + ((ischemia) * (GNa_GCaL_ICZ-GNa_GCaL_NZ)));
			}
	///////////////////////////////////////

real GNa_P = GNa * 1.25; // 11.7802*1.7; // BetaAdrenergic;

real fINap = camk_f_RyR;
real fINaLp = camk_f_RyR;
real fINa_P = fINaP; // PKA-P fraction as assigned as input, take the value 0 or 1
real fINa_BP = fINap*fINa_P;
real fINa_CaMKonly = fINap-fINa_BP;
real fINa_PKAonly = fINa_P-fINa_BP;

real INaBase_NP =  GNa*pow(m,3.0)*h*j; // Non-Phosphorylated 
real INaBase_CaMK = GNa*pow(m,3.0)*hp*jp;
real INaBase_PKA = GNa_P*pow(m_P,3.0)*h_P*j_P;
real INaBase_BP = GNa_P*pow(m_P,3.0)*hp_P*jp_P;

// 4 population 
real INaBase =  ((1-fINa_CaMKonly-fINa_PKAonly-fINa_BP)*INaBase_NP + fINa_CaMKonly*INaBase_CaMK + fINa_PKAonly*INaBase_PKA + fINa_BP*INaBase_BP) ;
real I_NaFast_junc = Fjunc * INaBase*(v-ena_junc);
real I_NaFast_sl = (1 - Fjunc) * INaBase * (v-ena_sl);

// INaL
//function [I_NaL_junc, I_NaL_sl, dmL, dhL, dhLp] = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ena_junc, ena_sl, layer, Fjunc, Fsl, INaL_Multiplier)
// calculate INaL

real mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
real tm = 0.1292 * exp(-pow((v+45.79)/15.54,2)) + 0.06487 * exp(-pow((v-4.823)/51.12,2));
real tmL=tm;
real dmL=(mLss-mL)/tmL;             // Rush-Larsen
real hLss=1.0/(1.0+exp((v+87.61)/7.488));
real thL=150.0;
real dhL=(hLss-hL)/thL;             // Rush-Larsen
real hLssp=1.0/(1.0+exp((v+93.81)/7.488));
real thLp=3.0*thL;
real dhLp=(hLssp-hLp)/thLp;
real GNaL=0.04229 * INaL_Multiplier * (1+fINaLp); // JAKUB added, based e.g. on https://physoc.onlinelibrary.wiley.com/doi/full/10.1113/EP085990
if (layer == EPI)
    GNaL=GNaL*0.6;

real I_NaL_junc = Fjunc * GNaL*(v-ena_junc)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
real I_NaL_sl = Fsl * GNaL*(v-ena_sl)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

// INa together aggregation variables:
// 1) for ionic update purposes, total sodium junc and sl are calculated
real I_Na_junc = I_NaFast_junc + I_NaL_junc;
real I_Na_sl = I_NaFast_sl + I_NaL_sl;

// 2) for plotting purposes, storing total fast and late sodium current
real I_NaFast = I_NaFast_junc+I_NaFast_sl;
real I_NaL = I_NaL_junc+I_NaL_sl;

// I_nabk: Na Background Current
real GNaB = 2* 0.297e-3 * INab_Multiplier;    // [mS/uF]
real I_Nabk_junc = Fjunc*GNaB*(v-ena_junc);
real I_Nabk_sl = Fsl*GNaB*(v-ena_sl);
real I_Nabk = I_Nabk_junc+I_Nabk_sl;

// I_nak: Na/K Pump Current
real IbarNaK = 2.10774 * INaK_Multiplier;     // [uA/uF]
real KmNaip = 11;         // [mM]
real KmNaip_PKA = 8.4615;
real KmKo = 1.5;         // [mM]
real fnak = 0.75 + (0.00375- ((140 - ions_na_o)/50) * 0.001)*v ; // Varying the slope mainly based on https://rupress.org/jgp/article-pdf/94/3/539/1814046/539.pdf
real I_nak_junc_noPKA = Fjunc*IbarNaK*fnak*ions_k_o /(1+pow(KmNaip/ions_na_junc,4)) /(ions_k_o+KmKo);
real I_nak_sl_noPKA = Fsl*IbarNaK*fnak*ions_k_o /(1+pow(KmNaip/ions_na_sl,4)) /(ions_k_o+KmKo);    
real I_nak_junc_PKA = Fjunc*IbarNaK*fnak*ions_k_o /(1+pow(KmNaip_PKA/ions_na_junc,4)) /(ions_k_o+KmKo);
real I_nak_sl_PKA = Fsl*IbarNaK*fnak*ions_k_o /(1+pow(KmNaip_PKA/ions_na_sl,4)) /(ions_k_o+KmKo);   

real I_nak_junc = (1 - fINaK_PKA) * I_nak_junc_noPKA + fINaK_PKA * I_nak_junc_PKA;
real I_nak_sl = (1 - fINaK_PKA) * I_nak_sl_noPKA + fINaK_PKA * I_nak_sl_PKA;
real I_nak = I_nak_junc+I_nak_sl;

// IKr_ToRORd_MM
// Variant based on Lu-Vandenberg
// Extracting state vector
int b = 0; // no channels blocked in via the mechanism of specific MM states
real vfrt = v*F/(R*T);

// transition rates
// from c0 to c1 in l-v model,
real alpha = 0.1161 * exp(0.2990 * vfrt);
// from c1 to c0 in l-v/
real beta =  0.2442 * exp(-1.604 * vfrt);

// from c1 to c2 in l-v/
real alpha1 = 1.25 * 0.1235 ;
// from c2 to c1 in l-v/
real beta1 =  0.1911;

// from c2 to o/           c1 to o
real alpha2 = 0.0578 * exp(0.9710 * vfrt); // from o to c2/
real beta2 = 0.349e-3* exp(-1.062 * vfrt); 

// from o to i
real alphai = 0.2533 * exp(0.5953 * vfrt); 
// from i to o
real betai = 0.04568 * exp(-0.8209 * vfrt); 

// from c2 to i (from c1 in orig)
real alphac2ToI = 0.52e-4 * exp(1.525 * vfrt); 
// from i to c2
// betaItoC2 = 0.85e-8 * exp(-1.842 * vfrt); %
real betaItoC2 = (beta2 * betai * alphac2ToI)/(alpha2 * alphai); 
// transitions themselves
// for reason of backward compatibility of naming of an older version of a
// MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

real dc0 = c1 * beta - c0 * alpha; // Euler
real dc1 = c0 * alpha + c2*beta1 - c1*(beta+alpha1); // Euler
real dc2 = c1 * alpha1 + o*beta2 + i*betaItoC2 - c2 * (beta1 + alpha2 + alphac2ToI); // Euler
real delta_o = c2 * alpha2 + i*betai - o*(beta2+alphai); // Euler
real di = c2*alphac2ToI + o*alphai - i*(betaItoC2 + betai); // Euler

real GKr = 0.04534 * sqrt(ions_k_o/5) * IKr_Multiplier;
if (layer == EPI)
    GKr=GKr*1.1; // changed compared to baseline
else if (layer==MID)
    GKr=GKr*0.7; // changed a bit

real IKr = GKr * o  * (v-ek);

// I_ks: Slowly Activating K Current
real pNaK = 0.01833;
real eks = (1/FoRT)*log((ions_k_o+pNaK*ions_na_o)/(ions_k_i+pNaK*ions_na_i));

// IKs_Bartos
real kPKA_Iks = fIKs_PKA;
real gks_factor_SA = 2.97002 * IKs_Multiplier *(pow(0.2, (2*apicobasal - 1)));

if (layer == MID)
    gks_factor_SA = 0.5 * gks_factor_SA;

real gks_factor = 0.01; 
real P_g_0 = gks_factor*(0.2+0.2*kPKA_Iks); 
real P_g_max = gks_factor*(0.8+7*kPKA_Iks); 
real P_vh_0 = -1-10*kPKA_Iks; 
real P_vh_max = -12-9*kPKA_Iks; 
real P_tau_0 = 26+9*kPKA_Iks;
real P_tau_max = 40+4*kPKA_Iks; 

real caks_junc = ions_ca_junc;
real caks_sl = ions_ca_sl; // normal simulation

real gks_junc = P_g_0 + (P_g_max-P_g_0)/(1 + pow(150e-6/caks_junc,1.3)); // Regulated by PKA
real gks_sl = P_g_0 + (P_g_max-P_g_0)/(1 + pow(150e-6/caks_sl,1.3)); // Regulated by PKA
real VsXs_Ca_junc = P_vh_0 + (P_vh_max-P_vh_0)/(1 + pow(350e-6/caks_junc,4)); // Regulated by PKA
real xsss_junc = 1/(1+exp(-(v-VsXs_Ca_junc)/25));
real VsTs_Ca_junc = P_tau_0 + (P_tau_max-P_tau_0)/(1 + pow(150e-6/caks_junc,3)); // Regulated by PKA
real tauxs_junc = 2*(50+(50+350*exp(-(pow(v+30,2))/4000))*1/(1+exp(-(v+VsTs_Ca_junc)/10)));
real VsXs_Ca_sl = P_vh_0 + (P_vh_max-P_vh_0)/(1 + pow(350e-6/caks_sl,4)); // Regulated by PKA
real xsss_sl = 1/(1+exp(-(v-VsXs_Ca_sl)/25));
real VsTs_Ca_sl = P_tau_0 + (P_tau_max-P_tau_0)/(1 + pow(150e-6/caks_sl,3)); // Regulated by PKA
real tauxs_sl = 2*(50+(50+350*exp(-(pow(v+30,2))/4000))*1/(1+exp(-(v+VsTs_Ca_sl)/10)));
real d_iks_xs_junc = (xsss_junc-iks_xs_junc)/tauxs_junc;    // Rush-Larsen
real d_iks_xs_sl = (xsss_sl-iks_xs_sl)/tauxs_sl;            // Rush-Larsen
//     eks = ek;

real I_ks_junc = Fjunc*gks_factor_SA*gks_junc*pow(iks_xs_junc,2)*(v-eks);
real I_ks_sl = Fsl*gks_factor_SA*gks_sl*pow(iks_xs_sl,2)*(v-eks);

real I_ks = I_ks_junc+I_ks_sl;

// IKb: background potassium current
real xkb = 1.0/(1.0+exp(-(v-10.8968)/(23.9871)));
real GKb = 0.010879 * IKb_Multiplier; 
if (layer==EPI)
    GKb=GKb*0.6;
real I_Kb = GKb*xkb*(v-ek);

// I_to: 
// Transient Outward K Current (slow and fast components)
real fItop = camk_f_RyR;

real GtoSlow, GtoFast;  
if (layer == EPI) {
    GtoSlow=0.02036 * Itos_Multiplier; //epi
    GtoFast=0.29856 * Itof_Multiplier; //epi0.88
} 
else if (layer == MID) {
    GtoSlow= 0.04632 * Itos_Multiplier; //mid
    GtoFast=0.14928 * Itof_Multiplier;  //mid; %epi - 0.875 scaling based on https://www.ahajournals.org/doi/pdf/10.1161/01.RES.72.3.671
}
else {
    GtoSlow=0.07210 * Itos_Multiplier; //endo
    GtoFast=0.01276 * Itof_Multiplier; //endo
}

real xtoss = 1/(1+exp(-(v-19.0)/13));
real ytoss = 1/(1+exp((v+19.5)/5));
real tauxtos = 9/(1+exp((v+3.0)/15))+0.5;
real tauytos = 800/(1+exp((v+60.0)/10))+30;
real d_ito_xtos = (xtoss - ito_xtos)/tauxtos;   // Rush-Larsen
real d_ito_ytos = (ytoss - ito_ytos)/tauytos;   // Rush-Larsen

real tauxtof = 8.5*exp(-pow((v+45)/50,2))+0.5;
real tauytof = 85*exp((-pow(v+40,2)/220))+7;
real d_ito_xtof = (xtoss - ito_xtof)/tauxtof;   // Rush-Larsen
real d_ito_ytof = (ytoss - ito_ytof)/tauytof;   // Rush-Larsen

// CaMKII effect on ss activation
real xtoss_p = 1/(1+exp(-(v-29.0)/13));
real d_ito_xtos_p = (xtoss_p - ito_xtos_p)/tauxtos; // Rush-Larsen
real d_ito_xtof_p = (xtoss_p - ito_xtof_p)/tauxtof; // Rush-Larsen

// And CaMKII effect on ss inactivation
real dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
real dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
real tauytos_p = tauytos * dti_develop * dti_recover;
real tauytof_p = tauytof * dti_develop * dti_recover;

real d_ito_ytos_p = (ytoss - ito_ytos_p)/tauytos_p; // Rush-Larsen
real d_ito_ytof_p = (ytoss - ito_ytof_p)/tauytof_p; // Rush-Larsen

real I_tos = GtoSlow*(v-ek) * ((1-fItop)*ito_xtos*ito_ytos + fItop*ito_xtos_p*ito_ytos_p);    // [uA/uF]
real I_tof = GtoFast*(v-ek)*((1-fItop) * ito_xtof*ito_ytof + fItop * ito_xtof_p*ito_ytof_p);
real I_to = I_tos + I_tof;

// I_k1: Time-Independent K Current
// IK1
real aK1 = 4.094/(1+exp(0.1217*(v-ek-49.934)));
real bK1 = (15.72*exp(0.0674*(v-ek-3.257))+exp(0.0618*(v-ek-594.31)))/(1+exp(-0.1629*(v-ek+14.207)));
real K1ss = aK1/(aK1+bK1);

real GK1 = IK1_Multiplier * 0.6992; 
if (layer==EPI)
    GK1=GK1*1.2;
else if (layer==MID)
    GK1=GK1*1.3;
real IK1 = GK1*sqrt(ions_k_o/5)*K1ss*(v-ek);

// I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
real GClCa = IClCa_Multiplier * 0.01615;   // [mS/uF]
real GClB = IClb_Multiplier * 0.00241;     // [mS/uF]
real KdClCa = 100e-3;    // [mM]

real I_ClCa_junc = 0.5 * Fjunc*GClCa/(1+KdClCa/ions_ca_junc)*(v-ecl);
real I_ClCa_sl = 0.5 * Fsl*GClCa/(1+KdClCa/ions_ca_sl)*(v-ecl);
real I_ClCa = I_ClCa_junc+I_ClCa_sl;
real I_Clbk = GClB*(v-ecl);

// I_Ca: L-type Calcium Current
real fICaLp = camk_f_ICaL;

real vffrt;
real Aff;
real Afs;
real f;
real Afcaf;
real Afcas;
real fca;
real fp;
real fcap;
real Kmn;
real k2n;
real km2n;
real anca;
real anca_i;
real Io;
real Ii;
real dielConstant;
real temp;
real constA;
real gamma_cai;
real gamma_cao;
real gamma_nai;
real gamma_nao;
real gamma_ki;
real gamma_ko;
real PhiCaL_ss;
real PhiCaNa_ss;
real PhiCaK_ss;
real gammaCaoMyo;
real gammaCaiMyo;
real PhiCaL_i;
real PhiCaNa_i;
real PhiCaK_i;
real sigmoidTransition, sigmoidTransition2;
real tauTransition, tauTransition2;
real PCap;
real PCaNa;
real PCaK;
real PCaNap;
real PCaKp;
real f_P;
real fcap_P;
real PCa_P;
real PCaNa_P;
real PCaK_P;
real fBP;
real fcaBP;
real fICaLP;
real fICaL_P;
real fICaL_BP;
real fICaL_CaMKonly;
real fICaL_PKAonly;
real ICaL_ss_NP;
real ICaL_ss_CaMK;
real ICaL_ss_PKA;
real ICaL_ss_BP;
real ICaL_i_NP;
real ICaL_i_CaMK;
real ICaL_i_PKA;
real ICaL_i_BP;
real ICaNa_ss_NP;
real ICaNa_ss_CaMK;
real ICaNa_ss_PKA;
real ICaNa_ss_BP;
real ICaNa_i_NP;
real ICaNa_i_CaMK;
real ICaNa_i_PKA;
real ICaNa_i_BP;
real ICaK_ss_NP;
real ICaK_ss_CaMK;
real ICaK_ss_PKA;
real ICaK_ss_BP;
real ICaK_i_NP;
real ICaK_i_CaMK;
real ICaK_i_PKA;
real ICaK_i_BP;
real PCa;
real rateRecovery;

real ICaL_junc; 
real ICaNa_junc;
real ICaK_junc;

real ICaL_sl;
real ICaNa_sl;
real ICaK_sl;

real dd;
real dff;
real dfs;
real dfcaf;
real dfcas;
real djca;
real dnca;
real dnca_i;
real dffp;
real dfcafp;
real d_ical_pureCDI_junc;
real d_ical_pureCDI_sl;
real dd_P;
real dff_P;
real dfs_P;
real dfcaf_P;
real dfcas_P;
real dfBPf;
real dfcaBPf;
real dss;
real fss;
real fcass;
real td;

real tff;
real tfs;
real tfcaf;
real tfcas;
real tjca;
real tffp;
real tfcafp;
real jcass;
real dPss;
real fss_P;
real fcass_P;
real fBPss;
real fcaBPss;


if (fICaL_PKA > 0) {
    vffrt=v*F*F/(R*T);

    //calculate ICaL, ICaNa, ICaK

    dss = fminf(1.0763*exp(-1.0070*exp(-0.0829*(v+3.62483))),1); 
    // if(v >31.4978) % activation cannot be greater than 1
    //     dss = 1;
    // end

    td = 1.5+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));

    dd = (dss-d)/td;       // Rush-Larsen
    fss = 1.0/(1.0+exp((v+19.58)/3.696));

    tff = 6.17111+1.0/(0.00126*exp(-(v+26.63596)/(9.69961))+0.00126*exp((v+26.63596)/(9.69961)));
    tfs = 2719.22489+1.0/(7.19411e-05*exp(-(v+5.74631)/(10.87690))+7.19411e-05*exp((v+5.74631)/(16.31535)));

    Aff=0.52477; 
    Afs=1.0-Aff;
    dff=(fss-ff)/tff;      // Rush-Larsen
    dfs=(fss-fs)/tfs;      // Rush-Larsen
    f=Aff*ff+Afs*fs;
    fcass=fss;
    tfcaf = 13.50673+1.0/(0.15420*exp(-(v-1.31611)/(11.33960))+0.15420*exp((v-1.31611)/(11.33960)));
    tfcas = 177.95813+1.0/(4.73955e-04*exp((-v+ 0.79049)/(0.81777)) + 4.73955e-04*exp((v+2.40474)/(1.90812)));

    Afcaf=0.3+0.6/(1.0+exp((v-9.24247)/(27.96201)));

    Afcas=1.0-Afcaf;
    dfcaf=(fcass-fcaf)/tfcaf;
    dfcas=(fcass-fcas)/tfcas;
    fca=Afcaf*fcaf+Afcas*fcas;

    tjca = 73.90298;
    jcass = 1.0/(1.0+exp((v+17.66945)/(3.21501)));
    djca = (jcass-jca)/tjca;   // Rush-Larsen
    tffp=2.5*tff;
    dffp=(fss-ffp)/tffp;
    fp=Aff*ffp+Afs*fs;
    tfcafp=2.5*tfcaf;
    dfcafp=(fcass-fcafp)/tfcafp;
    fcap=Afcaf*fcafp+Afcas*fcas;

    // SS nca
    Kmn=0.00222;
    k2n=957.85903;
    km2n=jca*0.84191;
    anca=1.0/(k2n/km2n+pow(1.0+Kmn/ions_ca_junc,3.80763));
    dnca=anca*k2n-nca*km2n;

    // myoplasmic nca
    anca_i = 1.0/(k2n/km2n+pow(1.0+Kmn/ions_ca_sl,3.80763));
    dnca_i = anca_i*k2n-nca_i*km2n;

    // SS driving force
    Io = 0.5*(ions_na_o + ions_k_o + ions_cl_o + 4*ions_ca_o)/1000; // ionic strength outside. /1000 is for things being in micromolar
    Ii = 0.5*(ions_na_sl + ions_k_i + ions_cl_i + 4*ions_ca_sl)/1000 ; // ionic strength inside. /1000 is for things being in micromolar % USING NONLOCAL CA
    // The ionic strength is too high for basic DebHuc. We'll use Davies
    dielConstant = 74; // water at 37°.
    temp = 310; // body temp in kelvins.
    constA = 1.82*pow(10,6)*pow(dielConstant*temp,-1.5);

    gamma_cai = pow(10,(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))); 
    gamma_cao = pow(10,(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_nai = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_nao = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_ki = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_ko = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));

    PhiCaL_ss =  4.0*vffrt*(gamma_cai*ions_ca_sl*exp(2.0*vfrt)-gamma_cao*ions_ca_o)/(exp(2.0*vfrt)-1.0); // USING NONLOCAL CA
    PhiCaNa_ss =  1.0*vffrt*(gamma_nai*ions_na_sl*exp(1.0*vfrt)-gamma_nao*ions_na_o)/(exp(1.0*vfrt)-1.0);
    PhiCaK_ss =  1.0*vffrt*(gamma_ki*ions_k_i*exp(1.0*vfrt)-gamma_ko*ions_k_o)/(exp(1.0*vfrt)-1.0);

    // sl/myo driving force
    Io = 0.5*(ions_na_o + ions_k_o + ions_cl_o + 4*ions_ca_o)/1000; // ionic strength outside. /1000 is for things being in micromolar
    Ii = 0.5*(ions_na_i + ions_k_i + ions_cl_i + 4*ions_ca_i)/1000; // ionic strength inside. /1000 is for things being in micromolar % USING NONLOCAL CA
    // The ionic strength is too high for basic DebHuc. We'll use Davies
    dielConstant = 74; // water at 37°.
    temp = 310; // body temp in kelvins.
    constA = 1.82*pow(10,6)*pow(dielConstant*temp,-1.5);

    gamma_cai = pow(10,(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))); 
    gamma_cao = pow(10,(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_nai = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_nao = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_ki = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_ko = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));

    gammaCaoMyo = gamma_cao;
    gammaCaiMyo = gamma_cai;

    PhiCaL_i =  4.0*vffrt*(gamma_cai*ions_ca_i*exp(2.0*vfrt)-gamma_cao*ions_ca_o)/(exp(2.0*vfrt)-1.0); // USING NONLOCAL CA
    PhiCaNa_i =  1.0*vffrt*(gamma_nai*ions_na_i*exp(1.0*vfrt)-gamma_nao*ions_na_o)/(exp(1.0*vfrt)-1.0);
    PhiCaK_i =  1.0*vffrt*(gamma_ki*ions_k_i*exp(1.0*vfrt)-gamma_ko*ions_k_o)/(exp(1.0*vfrt)-1.0);

    // Calculating "pure" CDI
    sigmoidTransition = 1 - 1/(1 + (1.86532 * ions_ca_junc/0.032));
    tauTransition = 1.09670 + (1-sigmoidTransition)*141.42990;

    rateRecovery = 0.02313;

    d_ical_pureCDI_junc = -ical_pureCDI_junc*sigmoidTransition/tauTransition + (1-ical_pureCDI_junc) * rateRecovery;

    // and for nondyadic
    sigmoidTransition2 = 1 - 1/(1 + (1.86532 * ions_ca_sl/0.032));
    tauTransition2 = 1.09670 + (1-sigmoidTransition2)*141.42990;

    d_ical_pureCDI_sl = -ical_pureCDI_sl*sigmoidTransition2/tauTransition2 + (1-ical_pureCDI_sl) * rateRecovery;
    // The rest
    real PCa = 1.5459e-04 * ICaL_PCaMultiplier;

    if (layer==EPI)
        PCa=PCa*1.025;
    else if (layer==MID)
        PCa=PCa*1.1;



	///////////////////   ICaL //////////////////////////////
	if(ischemia == 0) 	// NZ
	{
		PCa = PCa * GNa_GCaL_NZ;
	} else
	
		if(ischemia >= 1) 	// ICZ
		{
			PCa = PCa * GNa_GCaL_ICZ;
		} else
			// BZ, gradient
			{
				PCa = PCa * (GNa_GCaL_NZ + (ischemia * (GNa_GCaL_ICZ-GNa_GCaL_NZ)));
			}
	/////////////////////////////////////////////////////////////////
	
    PCap=1.1*PCa;
    PCaNa=1.1737/1.8969*0.00125*PCa;
    PCaK=1.1737/1.8969*3.574e-4*PCa;
    PCaNap=1.1737/1.8969*0.00125*PCap;
    PCaKp=1.1737/1.8969*3.574e-4*PCap;

    // PKA phosphorylation and other ones
    dPss=1.0323*exp(-1.0553*exp(-0.0810*(v+12.62483)));     // BetaAdrenergic ,fitted curve
    dPss=fminf(dPss,1);
    dd_P =(dPss-d_P)/td;    // Rush-Larsen
    fss_P=1.0/(1.0+exp((v+19.58+6)/3.696)); //BetaAdrenergic; 
    dff_P=(fss_P-ff_P)/tff;       // Rush-Larsen
    dfs_P=(fss_P-fs_P)/tfs;       // Rush-Larsen
    fcass_P=fss_P;
    dfcaf_P=(fcass_P-fcaf_P)/tfcaf;    // Rush-Larsen
    dfcas_P=(fcass_P-fcas_P)/tfcas;    // Rush-Larsen
    f_P=Aff*ff_P+Afs*fs_P;
    fcap_P=Afcaf*fcaf_P+Afcas*fcas_P;

    // The rest (PKA)
    PCa_P = PCa * 1.9;  //BetaAdrenergic;  (%Gong--> 1.9)
    if (layer==EPI)
        PCa_P=PCa_P*1.025;
    else if (layer==MID)
        PCa_P=PCa_P*1.1;

    PCaNa_P=1.1737/1.8969*0.00125*PCa_P;
    PCaK_P=1.1737/1.8969*3.574e-4*PCa_P;

    // Both-P population takes on dP, for f gate, takes ss from PKA and tau from CaMK (only fast component was modified)
    fBPss = fss_P;
    dfBPf = (fBPss-fBPf)/tffp;     // Rush-Larsen
    fBP = Aff*fBPf+Afs*fs_P; // only fast component modified, slow component same as fPs
    fcaBPss = fcass_P;
    dfcaBPf = (fcaBPss-fcaBPf)/tfcafp;  // Rush-Larsen
    fcaBP = Afcaf*fcaBPf+Afcas*fcas_P;
    
    // Get the PKA-P parameter value for this current
    fICaLP = fICaL_PKA;
    // fICaLp=fICaLp; % CaMK-P fraction
    fICaL_P = fICaLP; // PKA-P fraction
    fICaL_BP = fICaLp*fICaL_P;
    fICaL_CaMKonly = fICaLp-fICaL_BP;
    fICaL_PKAonly = fICaL_P-fICaL_BP;

    // ---------------------------------------------------

    ICaL_ss_NP=PCa*PhiCaL_ss*d*(f*(1.0-nca)+jca*fca*nca);
    ICaL_ss_CaMK=PCap*PhiCaL_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
    ICaL_ss_PKA=PCa_P*PhiCaL_ss*d_P*(f_P*(1.0-nca)+jca*fcap_P*nca);
    ICaL_ss_BP=PCa_P*PhiCaL_ss*d_P*(fBP*(1.0-nca)+jca*fcaBP*nca);
    ICaL_i_NP=PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
    ICaL_i_CaMK=PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
    ICaL_i_PKA=PCa_P*PhiCaL_i*d_P*(f_P*(1.0-nca_i)+jca*fcap_P*nca_i);
    ICaL_i_BP=PCa_P*PhiCaL_i*d_P*(fBP*(1.0-nca_i)+jca*fcaBP*nca_i);
    
    ICaNa_ss_NP=PCaNa*PhiCaNa_ss*d*(f*(1.0-nca)+jca*fca*nca);
    ICaNa_ss_CaMK=PCaNap*PhiCaNa_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
    ICaNa_ss_PKA=PCaNa_P*PhiCaNa_ss*d_P*(f_P*(1.0-nca)+jca*fcap_P*nca);
    ICaNa_ss_BP=PCaNa_P*PhiCaNa_ss*d_P*(fBP*(1.0-nca)+jca*fcaBP*nca);
    ICaNa_i_NP=PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
    ICaNa_i_CaMK=PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
    ICaNa_i_PKA=PCaNa_P*PhiCaNa_i*d_P*(f_P*(1.0-nca_i)+jca*fcap_P*nca_i);
    ICaNa_i_BP=PCaNa_P*PhiCaNa_i*d_P*(fBP*(1.0-nca_i)+jca*fcaBP*nca_i);

    ICaK_ss_NP=PCaK*PhiCaK_ss*d*(f*(1.0-nca)+jca*fca*nca);
    ICaK_ss_CaMK=PCaKp*PhiCaK_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
    ICaK_ss_PKA=PCaK_P*PhiCaK_ss*d_P*(f_P*(1.0-nca)+jca*fcap_P*nca);
    ICaK_ss_BP=PCaK_P*PhiCaK_ss*d_P*(fBP*(1.0-nca)+jca*fcaBP*nca);
    ICaK_i_NP=PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
    ICaK_i_CaMK=PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
    ICaK_i_PKA=PCaK_P*PhiCaK_i*d_P*(f_P*(1.0-nca_i)+jca*fcap_P*nca_i);
    ICaK_i_BP=PCaK_P*PhiCaK_i*d_P*(fBP*(1.0-nca_i)+jca*fcaBP*nca_i);

    // 4 population combination%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% And adding the pure
    // CDI modulation
    ICaL_junc = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaL_ss_NP + fICaL_CaMKonly*ICaL_ss_CaMK + fICaL_PKAonly*ICaL_ss_PKA + fICaL_BP*ICaL_ss_BP) * ical_pureCDI_junc;
    ICaNa_junc = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaNa_ss_NP + fICaL_CaMKonly*ICaNa_ss_CaMK + fICaL_PKAonly*ICaNa_ss_PKA + fICaL_BP*ICaNa_ss_BP) * ical_pureCDI_junc;
    ICaK_junc = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaK_ss_NP + fICaL_CaMKonly*ICaK_ss_CaMK + fICaL_PKAonly*ICaK_ss_PKA + fICaL_BP*ICaK_ss_BP) * ical_pureCDI_junc;
    
    ICaL_sl = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaL_i_NP + fICaL_CaMKonly*ICaL_i_CaMK + fICaL_PKAonly*ICaL_i_PKA + fICaL_BP*ICaL_i_BP) * ical_pureCDI_sl;
    ICaNa_sl = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaNa_i_NP + fICaL_CaMKonly*ICaNa_i_CaMK + fICaL_PKAonly*ICaNa_i_PKA + fICaL_BP*ICaNa_i_BP) * ical_pureCDI_sl;
    ICaK_sl = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaK_i_NP + fICaL_CaMKonly*ICaK_i_CaMK + fICaL_PKAonly*ICaK_i_PKA + fICaL_BP*ICaK_i_BP) * ical_pureCDI_sl;

    // And we weight ICaL (in ss) and ICaL_i
    ICaL_sl = ICaL_sl * (1-ICaL_fractionSS);
    ICaNa_sl = ICaNa_sl * (1-ICaL_fractionSS);
    ICaK_sl = ICaK_sl * (1-ICaL_fractionSS);
    ICaL_junc = ICaL_junc * ICaL_fractionSS;
    ICaNa_junc = ICaNa_junc * ICaL_fractionSS;
    ICaK_junc = ICaK_junc * ICaL_fractionSS;
}
else {
    vffrt=v*F*F/(R*T);

    // calculate ICaL, ICaNa, ICaK
    dss=fminf(1.0763*exp(-1.0070*exp(-0.0829*(v+3.62483))), 1); 
    // if(v >31.4978) % activation cannot be greater than 1
    //     dss = 1;
    // end

    td= 1.5+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));

    dd = (dss-d)/td; // Rush-Larsen
    fss = 1.0/(1.0+exp((v+19.58)/3.696));

    tff = 6.17111+1.0/(0.00126*exp(-(v+26.63596)/(9.69961))+0.00126*exp((v+26.63596)/(9.69961)));
    tfs = 2719.22489+1.0/(7.19411e-05*exp(-(v+5.74631)/(10.87690))+7.19411e-05*exp((v+5.74631)/(16.31535)));

    Aff=0.52477; 
    Afs=1.0-Aff;
    dff=(fss-ff)/tff;               // Rush-Larsen
    dfs=(fss-fs)/tfs;               // Rush-Larsen
    f=Aff*ff+Afs*fs;
    fcass=fss;
    tfcaf=13.50673+1.0/(0.15420*exp(-(v-1.31611)/(11.33960))+0.15420*exp((v-1.31611)/(11.33960)));
    tfcas=177.95813+1.0/(4.73955e-04*exp((-v+ 0.79049)/(0.81777)) + 4.73955e-04*exp((v+2.40474)/(1.90812)));

    Afcaf = 0.3+0.6/(1.0+exp((v-9.24247)/(27.96201)));

    Afcas=1.0-Afcaf;
    dfcaf=(fcass-fcaf)/tfcaf;       // Rush-Larsen
    dfcas=(fcass-fcas)/tfcas;       // Rush-Larsen
    fca=Afcaf*fcaf+Afcas*fcas;

    tjca = 73.90298;
    jcass = 1.0/(1.0+exp((v+17.66945)/(3.21501)));
    djca=(jcass-jca)/tjca;          // Rush-Larsen
    tffp=2.5*tff;
    dffp=(fss-ffp)/tffp;            // Rush-Larsen
    fp=Aff*ffp+Afs*fs;
    tfcafp=2.5*tfcaf;
    dfcafp=(fcass-fcafp)/tfcafp;    // Rush-Larsen
    fcap=Afcaf*fcafp+Afcas*fcas;

    // SS nca
    Kmn=0.00222;
    k2n=957.85903;
    km2n=jca*0.84191;
    anca=1.0/(k2n/km2n+pow(1.0+Kmn/ions_ca_junc,3.80763));
    dnca=anca*k2n-nca*km2n;

    // myoplasmic nca
    anca_i = 1.0/(k2n/km2n+pow(1.0+Kmn/ions_ca_sl,3.80763));
    dnca_i = anca_i*k2n-nca_i*km2n;

    // SS driving force
    Io = 0.5*(ions_na_o + ions_k_o + ions_cl_o + 4*ions_ca_o)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
    Ii = 0.5*(ions_na_sl + ions_k_i + ions_cl_i + 4*ions_ca_sl)/1000 ; // ionic strength inside. /1000 is for things being in micromolar % USING NONLOCAL CA
    // The ionic strength is too high for basic DebHuc. We'll use Davies
    dielConstant = 74; // water at 37°.
    temp = 310; // body temp in kelvins.
    constA = 1.82*pow(10,6)*pow(dielConstant*temp,-1.5);

    gamma_cai = pow(10,(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))); 
    gamma_cao = pow(10,(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_nai = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_nao = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_ki = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_ko = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));

    PhiCaL_ss =  4.0*vffrt*(gamma_cai*ions_ca_sl*exp(2.0*vfrt)-gamma_cao*ions_ca_o)/(exp(2.0*vfrt)-1.0); // USING NONLOCAL CA
    PhiCaNa_ss =  1.0*vffrt*(gamma_nai*ions_na_sl*exp(1.0*vfrt)-gamma_nao*ions_na_o)/(exp(1.0*vfrt)-1.0);
    PhiCaK_ss =  1.0*vffrt*(gamma_ki*ions_k_i*exp(1.0*vfrt)-gamma_ko*ions_k_o)/(exp(1.0*vfrt)-1.0);

    // sl/myo driving force
    Io = 0.5*(ions_na_o + ions_k_o + ions_cl_o + 4*ions_ca_o)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
    Ii = 0.5*(ions_na_i + ions_k_i + ions_cl_i + 4*ions_ca_i)/1000 ; // ionic strength inside. /1000 is for things being in micromolar % USING NONLOCAL CA
    // The ionic strength is too high for basic DebHuc. We'll use Davies
    dielConstant = 74; // water at 37°.
    temp = 310; // body temp in kelvins.
    constA = 1.82*pow(10,6)*pow(dielConstant*temp,-1.5);

    gamma_cai = pow(10,(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))); 
    gamma_cao = pow(10,(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_nai = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_nao = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));
    gamma_ki = pow(10,(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii)));
    gamma_ko = pow(10,(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io)));

    gammaCaoMyo = gamma_cao;
    gammaCaiMyo = gamma_cai;

    PhiCaL_i =  4.0*vffrt*(gamma_cai*ions_ca_i*exp(2.0*vfrt)-gamma_cao*ions_ca_o)/(exp(2.0*vfrt)-1.0); // USING NONLOCAL CA
    PhiCaNa_i =  1.0*vffrt*(gamma_nai*ions_na_i*exp(1.0*vfrt)-gamma_nao*ions_na_o)/(exp(1.0*vfrt)-1.0);
    PhiCaK_i =  1.0*vffrt*(gamma_ki*ions_k_i*exp(1.0*vfrt)-gamma_ko*ions_k_o)/(exp(1.0*vfrt)-1.0);

    // Calculating "pure" CDI
    sigmoidTransition = 1 - 1/(1 + pow(1.86532 * ions_ca_junc/0.032,0.99613));
    tauTransition = 1.09670 + (1-sigmoidTransition)*141.42990;

    rateRecovery = 0.02313;

    d_ical_pureCDI_junc = -ical_pureCDI_junc*sigmoidTransition/tauTransition + (1-ical_pureCDI_junc) * rateRecovery;

    // and for nondyadic
    sigmoidTransition2 = 1 - 1/(1 + pow(1.86532 * ions_ca_sl/0.032,0.99613));
    tauTransition2 = 1.09670 + (1-sigmoidTransition2)*141.42990;

    d_ical_pureCDI_sl = -ical_pureCDI_sl*sigmoidTransition2/tauTransition2 + (1-ical_pureCDI_sl) * rateRecovery;
    // The rest
    PCa=1.5459e-04 * ICaL_PCaMultiplier;
    
    if (layer==EPI)
        PCa=PCa*1.025;
    else if (layer==MID)
        PCa=PCa*1.1;

    PCap=1.1*PCa;
    PCaNa=1.1737/1.8969*0.00125*PCa;
    PCaK=1.1737/1.8969*3.574e-4*PCa;
    PCaNap=1.1737/1.8969*0.00125*PCap;
    PCaKp=1.1737/1.8969*3.574e-4*PCap;

    dd_P=0;
    dff_P=0;
    dfs_P=0;
    dfcaf_P=0;      //(fcass_P-fcaf_P)/tfcaf;
    dfcas_P=0;      //(fcass_P-fcas_P)/tfcas;
    dfBPf = 0;      //(fBPss-fBPf)/tffp;
    dfcaBPf = 0;    // (fcaBPss-fcaBPf)/tfcafp;  
    fICaL_CaMKonly = fICaLp;

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ICaL_ss_NP=PCa*PhiCaL_ss*d*(f*(1.0-nca)+jca*fca*nca);
    ICaL_ss_CaMK=PCap*PhiCaL_ss*d*(fp*(1.0-nca)+jca*fcap*nca);

    ICaL_i_NP=PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
    ICaL_i_CaMK=PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
    
    ICaNa_ss_NP=PCaNa*PhiCaNa_ss*d*(f*(1.0-nca)+jca*fca*nca);
    ICaNa_ss_CaMK=PCaNap*PhiCaNa_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
    ICaNa_i_NP=PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
    ICaNa_i_CaMK=PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);

    real ICaK_ss_NP=PCaK*PhiCaK_ss*d*(f*(1.0-nca)+jca*fca*nca);
    real ICaK_ss_CaMK=PCaKp*PhiCaK_ss*d*(fp*(1.0-nca)+jca*fcap*nca);

    real ICaK_i_NP=PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
    real ICaK_i_CaMK=PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);

    // 4 population combination 
    // And adding the pure
    // CDI modulation
    ICaL_junc = ((1-fICaL_CaMKonly)*ICaL_ss_NP + fICaL_CaMKonly*ICaL_ss_CaMK) * ical_pureCDI_junc;
    ICaNa_junc = ((1-fICaL_CaMKonly)*ICaNa_ss_NP + fICaL_CaMKonly*ICaNa_ss_CaMK ) * ical_pureCDI_junc;
    ICaK_junc = ((1-fICaL_CaMKonly)*ICaK_ss_NP + fICaL_CaMKonly*ICaK_ss_CaMK) * ical_pureCDI_junc;
    
    ICaL_sl = ((1-fICaL_CaMKonly)*ICaL_i_NP + fICaL_CaMKonly*ICaL_i_CaMK) * ical_pureCDI_sl;
    ICaNa_sl = ((1-fICaL_CaMKonly)*ICaNa_i_NP + fICaL_CaMKonly*ICaNa_i_CaMK) * ical_pureCDI_sl;
    ICaK_sl = ((1-fICaL_CaMKonly)*ICaK_i_NP + fICaL_CaMKonly*ICaK_i_CaMK) * ical_pureCDI_sl;

    // And we weight ICaL (in ss) and ICaL_i
    ICaL_sl = ICaL_sl * (1-ICaL_fractionSS);
    ICaNa_sl = ICaNa_sl * (1-ICaL_fractionSS);
    ICaK_sl = ICaK_sl * (1-ICaL_fractionSS);
    ICaL_junc = ICaL_junc * ICaL_fractionSS;
    ICaNa_junc = ICaNa_junc * ICaL_fractionSS;
    ICaK_junc = ICaK_junc * ICaL_fractionSS;
}

real I_CaK = ICaK_junc + ICaK_sl;
real I_Ca = ICaL_junc+ICaL_sl;
real I_CaNa = ICaNa_junc+ICaNa_sl;
real I_Catot = I_Ca+I_CaK+I_CaNa;

// I_ncx: Na/Ca Exchanger flux
// function for sl and outer dyad ncx
//[ I_ncx_sl, I_ncx_junc] = getINaCa_ORd2011(membrane_v,F,R,T, ions_na_junc, ions_na_sl, ions_na_i, ions_na_o, ions_ca_junc, ions_ca_sl, ions_ca_i, ions_ca_o, layer, INaCa_Multiplier, INaCa_fractionSS);
real zca=2.0;
real kna1=15.0;
real kna2=5.0;
real kna3=88.12; 
real kasymm=12.5;
real wna=6.0e4;
real wca=6.0e4;
real wnaca=5.0e3;
real kcaon=1.5e6 * 1.25;
real kcaoff=5.0e3;
real qna=0.5224;
real qca=0.1670 ;
real hca=exp((qca*v*F)/(R*T));
real hna=exp((qna*v*F)/(R*T));
real h1=1+ions_na_i/kna3*(1+hna);
real h2=(ions_na_i*hna)/(kna3*h1);
real h3=1.0/h1;
real h4=1.0+ions_na_i/kna1*(1+ions_na_i/kna2);
real h5=ions_na_i*ions_na_i/(h4*kna1*kna2);
real h6=1.0/h4;
real h7=1.0+ions_na_o/kna3*(1.0+1.0/hna);
real h8=ions_na_o/(kna3*hna*h7);
real h9=1.0/h7;
real h10=kasymm+1.0+ions_na_o/kna1*(1.0+ions_na_o/kna2);
real h11=ions_na_o*ions_na_o/(h10*kna1*kna2);
real h12=1.0/h10;
real k1=h12*ions_ca_o*kcaon;
real k2=kcaoff;
real k3p=h9*wca;
real k3pp=h8*wnaca;
real k3=k3p+k3pp;
real k4p=h3*wca/hca;
real k4pp=h2*wnaca;
real k4=k4p+k4pp;
real k5=kcaoff;
real k6=h6*ions_ca_i*kcaon;
real k7=h5*h2*wna;
real k8=h8*h11*wna;
real x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
real x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
real x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
real x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
real E1=x1/(x1+x2+x3+x4);
real E2=x2/(x1+x2+x3+x4);
real E3=x3/(x1+x2+x3+x4);
real E4=x4/(x1+x2+x3+x4);
real KmCaAct=150e-6;
real allo=1.0/(1.0+pow(KmCaAct/ions_ca_sl,2));
// allo = 1 - 1/(1+((ions_ca_sl*1e6)/107)^2.57);
real zna=1.0;
real JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
real JncxCa=E2*k2-E1*k1;
real Gncx= 0.00179* INaCa_Multiplier;
if (layer==EPI)
    Gncx=Gncx*1.1;
else if (layer==MID)
    Gncx=Gncx*1.4;
real I_ncx_sl=(1-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa);

// calculate INaCa_ss
h1=1+ions_na_sl/kna3*(1+hna);
h2=(ions_na_sl*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+ions_na_sl/kna1*(1+ions_na_sl/kna2);
h5=ions_na_sl*ions_na_sl/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+ions_na_o/kna3*(1.0+1.0/hna);
h8=ions_na_o/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+ions_na_o/kna1*(1+ions_na_o/kna2);
h11=ions_na_o*ions_na_o/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*ions_ca_o*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*ions_ca_sl*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150e-6;
allo=1.0/(1.0+pow(KmCaAct/ions_ca_junc,2));
// allo = 1 - 1/(1+((ions_ca_junc*1e6)/107)^2.57);
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
real I_ncx_junc=INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa);

real I_ncx = I_ncx_junc+I_ncx_sl; 

// I_pca: Sarcolemmal Ca Pump Current
real IbarSLCaP = 0.02064*IpCa_Multiplier; 
real KmPCa = 0.5e-3;     // [mM]

real Q10SLCaP = 2.35;    // [none]
real I_pca_junc = Fjunc*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(ions_ca_junc,1.6)/(pow(KmPCa,1.6)+pow(ions_ca_junc,1.6));
real I_pca_sl = Fsl*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(ions_ca_sl,1.6)/(pow(KmPCa,1.6)+pow(ions_ca_sl,1.6));
real I_pca = I_pca_junc+I_pca_sl;

// I_cabk: Ca Background Current
real GCaB = ICab_Multiplier * 5.15575e-04;    // [uA/uF]
real I_cabk_junc = Fjunc*GCaB*(v-eca_junc);
real I_cabk_sl = Fsl*GCaB*(v-eca_sl);
real I_cabk = I_cabk_junc+I_cabk_sl;


 // IKATP
real A_atp = 2.0;
real K_atp = 0.25;
real K_o_n = 5.0;
real fkatp = 0.0;

//////////////////// I_Katp ////////////////////////////////////////////////////


if(ischemia == 0) 	// NZ
{
	fkatp = fKatp_NZ;
} else

	if(ischemia >= 1) 	// ICZ
	{
		fkatp = fKatp_ICZ;
	} else
		// BZ, gradient
		{
			fkatp = fKatp_NZ + (ischemia * (fKatp_ICZ-fKatp_NZ));
		}

real gkatp = 4.3195;
real akik = pow((ions_k_o/K_o_n),0.24);
real bkik = 1.0/(1.0 + pow((A_atp/K_atp),2));
real IKatp = fkatp * gkatp * akik * bkik * (v-ek);



// SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
//[J_SRCarel, J_SRleak, Jrel_ICaLdep, d_ryr_R, d_ryr_O, d_ryr_I, d_ryr_CaRI, d_ryr_R_p, d_ryr_O_p, d_ryr_I_p, d_ryr_CaRI_p, d_jrel_icaldep_act, d_jrel_icaldep_f1, d_jrel_icaldep_f2] = getJrel(ions_ca_junc, ions_ca_SR, I_Ca_junc, ryr_R, ryr_O, ryr_I, ryr_CaRI, ryr_R_p, ryr_O_p, ryr_I_p, jrel_icaldep_act, jrel_icaldep_f1, jrel_icaldep_f2, ryr_CaRI_p, camk_f_RyR, Jrel_Multiplier, extraParams);
real ks = 26.6 * Jrel_Multiplier;                 // [1/ms]
real koCa = 23.87221;               // [mM^-2 1/ms]   %default 10   modified 20
real kom = 0.16219;              // [1/ms]
real kiCa = 0.39871;              // [1/mM/ms]
real kim = 0.04311;             // [1/ms]
real ec50SR = 0.75385;           // [mM]
real steepnessCaSR = 5.09473;
real caExpFactor = 0.68655 ;
real caTransFactor = 0.94428 ;
real caExpFactor2 = 2.06273 ;
real caTransFactor2 = 0.52967 ;

// Directly ICaL-coupled RyRs
real directRelMidpoint = 0.95271;
real bt=12.47670; 
real a_rel=1.25 * bt ;
real I_Ca_junc_positive = abs(ICaL_junc);
real I_Ca_junc_sigmoided = 1 - 1/(1+pow(I_Ca_junc_positive/0.45,4.5));
real Jrel_inf=a_rel*(I_Ca_junc_sigmoided)/(1.0+pow(directRelMidpoint/ions_ca_SR,7.72672)); 
real tau_rel=bt/(1.0+0.0123/ions_ca_SR);

if (tau_rel<0.001)
    tau_rel=0.001;

real d_jrel_icaldep_act = (Jrel_inf-jrel_icaldep_act)/tau_rel;  // Rush-Larsen

// and inactivation
real tauInact = 64.11202;
real Jrel_inact_inf = 1/(1 + (I_Ca_junc_sigmoided/1e-3));
real d_jrel_icaldep_f1 = (Jrel_inact_inf - jrel_icaldep_f1)/tauInact;   // Rush-Larsen

// and slower inactivation
real tauInact2 = 119.48978;
real Jrel_inact_inf2 = 1/(1 + (I_Ca_junc_sigmoided/0.6e-3));
real d_jrel_icaldep_f2 = (Jrel_inact_inf2 - jrel_icaldep_f2)/tauInact2; // Rush-Larsen

real Jrel_ICaLdep = (0.00174 * jrel_icaldep_act * jrel_icaldep_f1 * jrel_icaldep_f2 ); 

// Main Ca-sensitive RyRs
real MaxSR = 15; 
real MinSR = 1;
real kCaSR = MaxSR - (MaxSR-MinSR)/(1+pow(ec50SR/ions_ca_SR,steepnessCaSR));
real koSRCa = koCa/kCaSR;
real kiSRCa = kiCa*kCaSR;

real ecCaI = 0.001;
real steepnessCaI = 5.93447; // change compared to 5.7.4.5.4 was reflected here
real minCaI = 0.93249;
real maxCaI = 30.13294;

real baseRateCaI = 3.02320e-04;
real sigmoidBaseCaI = minCaI + (maxCaI-minCaI)/(1+pow(ecCaI/ions_ca_junc,steepnessCaI));
real RI_to_CI = baseRateCaI * sigmoidBaseCaI;
real CI_to_RI = 0.00248; 

// not phosphorylated by CaMKII
real RIcleft = 1 - ryr_R - ryr_O - ryr_I - ryr_CaRI;
real d_ryr_R = (kim*RIcleft-kiSRCa*(caTransFactor*pow(ions_ca_junc,caExpFactor))*ryr_R)-(caTransFactor2*koSRCa*pow(ions_ca_junc,caExpFactor2)*ryr_R-kom*ryr_O);     // R
real d_ryr_O = (caTransFactor2*koSRCa*pow(ions_ca_junc,caExpFactor2)*ryr_R-kom*ryr_O)-(kiSRCa*(caTransFactor*pow(ions_ca_junc,caExpFactor))*ryr_O-kim*ryr_I);       // O
real d_ryr_I = (kiSRCa*(caTransFactor*pow(ions_ca_junc,caExpFactor))*ryr_O-kim*ryr_I)-(kom*ryr_I-caTransFactor2*koSRCa*pow(ions_ca_junc,caExpFactor2)*RIcleft);               // I
real d_ryr_CaRI =  RI_to_CI*RIcleft - CI_to_RI * ryr_CaRI; // shift of 29 state numbers

real J_SRCarel_np = ks*ryr_O*(ions_ca_SR-ions_ca_junc) + Jrel_ICaLdep;          // [mM/ms]

// And also a version of phosphorylated
caTransFactor2 = caTransFactor2 * 1.5; // *1.69
real RIcleftP = 1 - ryr_R_p - ryr_O_p - ryr_I_p - ryr_CaRI_p;
real d_ryr_R_p = (kim*RIcleftP-kiSRCa*(caTransFactor*pow(ions_ca_junc,caExpFactor))*ryr_R_p)-(caTransFactor2*koSRCa*pow(ions_ca_junc,caExpFactor2)*ryr_R_p-kom*ryr_O_p);   // R
real d_ryr_O_p = (caTransFactor2*koSRCa*pow(ions_ca_junc,caExpFactor2)*ryr_R_p-kom*ryr_O_p)-(kiSRCa*(caTransFactor*pow(ions_ca_junc,caExpFactor))*ryr_O_p-kim*ryr_I_p);     // O
real d_ryr_I_p = (kiSRCa*(caTransFactor*pow(ions_ca_junc,caExpFactor))*ryr_O_p-kim*ryr_I_p)-(kom*ryr_I_p-caTransFactor2*koSRCa*pow(ions_ca_junc,caExpFactor2)*RIcleftP);   // I
real d_ryr_CaRI_p =  RI_to_CI*RIcleftP - CI_to_RI * ryr_CaRI_p;

real J_SRCarel_p = ks*ryr_O_p*(ions_ca_SR-ions_ca_junc) + Jrel_ICaLdep;          // [mM/ms]

// Total release
real J_SRCarel = J_SRCarel_p * camk_f_RyR + J_SRCarel_np * (1 - camk_f_RyR);

// Additional leak
real nonlinearModifier = 0.2144*exp(1.83*ions_ca_SR);  // Leak should be nonlinear with load
real CaMKIILeakMultiplier = 1 + 2*camk_f_RyR;   // And it should be promoted by CaMKII
real J_SRleak = 1.59306e-06*(ions_ca_SR-ions_ca_junc) * nonlinearModifier * CaMKIILeakMultiplier;           //   [mM/ms]

// Reuptake SERCA
//J_serca = getJup(ions_ca_SR, ions_ca_i, Temp, camk_f_PLB, fPLB_PKA, casig_SERCA_act, layer, Jup_Multiplier, extraParams);
real Q10SRCaP = 2.6;          // [none]
real Vmax_SRCaP = 0.00543 * Jup_Multiplier;  // [mM/msec] (286 umol/L cytosol/sec)

real Kmr = 2.31442;               // [mM]L cytosol
real hillSRCaP =  1.02809;       // [mM]
real Kmf = 0.30672e-03;          // [mM] default

if (layer == EPI)
    Vmax_SRCaP = Vmax_SRCaP * 1.2; // based loosely on https://www.ahajournals.org/doi/10.1161/01.RES.62468.25308.27?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed - there may be other datasets, and other values. % the 20% is just an initial guesstimate.   
    
// PLB phosphorylation effect on affinity
real phosphorylationTotal = camk_f_PLB + fPLB_PKA - camk_f_PLB*fPLB_PKA; // we assume the same effect, just making sure we don't count it twice.
real Kmf_Phospho = Kmf * 0.5; //; % Similar percentage effect as in Heijman 2011 % CHANGED JAKUB

// Direct Ca-based acceleration
real Km_SERCA_Ca = 0.4; //0.03 in Heijman 2011; affinity for direct Vmax modulation by CaMKII
real Max_Vmax_SERCA_Ca = 1.11142; 
real Vmax_mult = 1 + Max_Vmax_SERCA_Ca / (1 + pow(Km_SERCA_Ca / casig_SERCA_act,2));
        
//     J_serca = Q10SRCaP^Qpow*Vmax_SRCaP*Vmax_mult*((ions_ca_i/Km_eff)^hillSRCaP-(ions_ca_SR/Kmr)^hillSRCaP)...
//         /(1+(ions_ca_i/Km_eff)^hillSRCaP+(ions_ca_SR/Kmr)^hillSRCaP);
real J_serca_np = pow(Q10SRCaP,Qpow)*Vmax_SRCaP*Vmax_mult*(pow((ions_ca_i/Kmf),hillSRCaP)-pow((ions_ca_SR/Kmr),hillSRCaP))/(1+pow((ions_ca_i/Kmf),hillSRCaP)+pow((ions_ca_SR/Kmr),hillSRCaP));
real J_serca_p = pow(Q10SRCaP,Qpow)*Vmax_SRCaP*Vmax_mult*(pow((ions_ca_i/Kmf_Phospho),hillSRCaP)-pow((ions_ca_SR/Kmr),hillSRCaP))/(1+pow((ions_ca_i/Kmf_Phospho),hillSRCaP)+pow((ions_ca_SR/Kmr),hillSRCaP));
real J_serca = J_serca_np * (1 - phosphorylationTotal) + J_serca_p * phosphorylationTotal;

// --------------------------------------------------------------------
// Land-Niederer model of contraction
//[Ta, d_contraction_XS, d_contraction_XW, d_contraction_Ca_TRPN, d_contraction_TmBlocked, d_contraction_ZETAS, d_contraction_ZETAW] = getContractionLand(ions_ca_i, contraction_XS, contraction_XW, contraction_Ca_TRPN, contraction_TmBlocked, contraction_ZETAS, contraction_ZETAW,  fTnI_PKA, fMyBPC_PKA, extraParams);
real fracTnIpo = .0031;  // Derived quantity (TnI_PKAp(baseline)/TnItot)
real fPKA_TnI = (1.45-0.45*(1-fTnI_PKA)/(1-fracTnIpo)); // multiplier for Ca unbinding from troponin.

real PKAForceMultiplier = 1 + fMyBPC_PKA * 0.26; 
real PKAXBacceleration = 1 + 0.5 * fMyBPC_PKA;

//EC parameters
real perm50=0.35;
real TRPN_n=1.65;
real koff= 0.07854;
real dr=0.25;
real wfrac=0.5;
real TOT_A=25;
real ktm_unblock= 0.02626; 
real beta_1=-2.4;
real beta_0=2.3;
real gamma=0.0085;
real gamma_wu=0.615;
real phi=2.23;

real nperm;
real ca50;
real Tref;
real nu;
real mu;
if (mode == 1) {
    nperm=2.2;
    ca50=2.5;
    Tref=40.5;
    nu=1;
    mu=1;
}
else {
    nperm=2.036;
    ca50=fPKA_TnI * 0.76450;
    Tref= 80; //120; % Jakub updated
    nu=10.15996;
    mu=3.94046;
}
 
real k_ws=0.004 * (1 + fMyBPC_PKA/2); 
real k_uw=0.026; 

real lambda_min=0.87;
real lambda_max=1.2;

k_ws=k_ws*mu;
k_uw=k_uw*nu;

real cdw=phi*k_uw*(1-dr)*(1-wfrac)/((1-dr)*wfrac);
real cds=phi*k_ws*(1-dr)*wfrac/dr;
real k_wu=k_uw*(1/wfrac-1)-k_ws;
real k_su=k_ws*(1/dr-1)*wfrac;
real A=(0.25*TOT_A)/((1-dr)*wfrac+dr)*(dr/0.25);

// XB model
real lambda0 = fminf(lambda_max,lambda);
real Lfac = fmaxf(0,1+beta_0*(lambda0+fminf(lambda_min,lambda0)-(1+lambda_min)));

real XU=(1-contraction_TmBlocked)-contraction_XW-contraction_XS; // unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap
real xb_ws=k_ws*contraction_XW ;
real xb_uw=k_uw*XU* (1 + fMyBPC_PKA/2);
real xb_wu=k_wu*contraction_XW;
real xb_su=k_su*contraction_XS;

real gamma_rate;
real term1 = (contraction_ZETAS > 0) ? contraction_ZETAS : 0.0f;
real term2 = (contraction_ZETAS < -1) ? (-contraction_ZETAS - 1) : 0.0f;
gamma_rate = gamma * fmaxf(term1, term2);
//gamma_rate=gamma*max((contraction_ZETAS>0).*contraction_ZETAS,(contraction_ZETAS<-1).*(-contraction_ZETAS-1));

real xb_su_gamma=gamma_rate*contraction_XS;
real gamma_rate_w = gamma_wu*abs(contraction_ZETAW); // weak xbs don't like being strained
real xb_wu_gamma=gamma_rate_w*contraction_XW;

real d_contraction_XS=xb_ws-xb_su-xb_su_gamma;
real d_contraction_XW=xb_uw-xb_wu-xb_ws-xb_wu_gamma;

ca50=ca50+beta_1*fminf(0.2,lambda-1);
real d_contraction_Ca_TRPN=koff*(pow((ions_ca_i*1000)/ca50,TRPN_n)*(1-contraction_Ca_TRPN) - contraction_Ca_TRPN); 

real XSSS=dr*0.5;
real XWSS=(1-dr)*wfrac*0.5;
real ktm_block=ktm_unblock*(pow(perm50,nperm))*0.5/(0.5-XSSS-XWSS);
real d_contraction_TmBlocked=ktm_block*fminf(100,(pow(contraction_Ca_TRPN,-(nperm/2))))*XU-ktm_unblock*(pow(contraction_Ca_TRPN,(nperm/2)))*contraction_TmBlocked;

//velocity dependence -- assumes distortion resets on W->S
real d_contraction_ZETAS=A*lambda_rate-cds*contraction_ZETAS;// - gamma_rate * ZETAS;
real d_contraction_ZETAW=A*lambda_rate-cdw*contraction_ZETAW;// - gamma_rate_w * ZETAW;

// Active Force
real Ta=Lfac*PKAForceMultiplier * (Tref/dr)*((contraction_ZETAS+1)*contraction_XS+(contraction_ZETAW)*contraction_XW);

// stimulus current
real Istim = stim_current;

// Sodium and Calcium Buffering
real d_buffers_NaBj = kon_na*ions_na_junc*(Bmax_Naj-buffers_NaBj)-koff_na*buffers_NaBj;         // NaBj      [mM/ms]
real d_buffers_NaBsl = kon_na*ions_na_sl*(Bmax_Nasl-buffers_NaBsl)-koff_na*buffers_NaBsl;       // NaBsl     [mM/ms]

real d_buffers_TnClow = Bmax_TnClow * d_contraction_Ca_TRPN; // myofilament - based on LAND contraction model
real d_buffers_TnCHc = kon_tnchca*ions_ca_i*(Bmax_TnChigh-buffers_TnCHc-buffers_TnCHm)-koff_tnchca*buffers_TnCHc; // TnCHc     [mM/ms]
real d_buffers_TnCHm = kon_tnchmg*ions_mg_i*(Bmax_TnChigh-buffers_TnCHc-buffers_TnCHm)-koff_tnchmg*buffers_TnCHm;   // TnCHm     [mM/ms]
real d_buffers_CaM = kon_cam*ions_ca_i*(Bmax_CaM-buffers_CaM)-koff_cam*buffers_CaM;                 // CaM       [mM/ms]
real d_buffers_Myosin_ca = kon_myoca*ions_ca_i*(Bmax_myosin-buffers_Myosin_ca-buffers_Myosin_mg)-koff_myoca*buffers_Myosin_ca;    // Myosin_ca [mM/ms]
real d_buffers_Myosin_mg = kon_myomg*ions_mg_i*(Bmax_myosin-buffers_Myosin_ca-buffers_Myosin_mg)-koff_myomg*buffers_Myosin_mg;    // Myosin_mg [mM/ms]
real d_buffers_SRB = kon_sr*ions_ca_i*(Bmax_SR-buffers_SRB)-koff_sr*buffers_SRB;                    // SRB       [mM/ms]

real J_CaB_cytosol = d_buffers_TnClow + d_buffers_TnCHc + d_buffers_CaM + d_buffers_Myosin_ca + d_buffers_SRB;

// Junctional and SL Ca Buffers
real d_buffers_SLLj = kon_sll*ions_ca_junc*(Bmax_SLlowj-buffers_SLLj)-koff_sll*buffers_SLLj;       // SLLj      [mM/ms]
real d_buffers_SLLsl = kon_sll*ions_ca_sl*(Bmax_SLlowsl-buffers_SLLsl)-koff_sll*buffers_SLLsl;      // SLLsl     [mM/ms]
real d_buffers_SLHj = kon_slh*ions_ca_junc*(Bmax_SLhighj-buffers_SLHj)-koff_slh*buffers_SLHj;      // SLHj      [mM/ms]
real d_buffers_SLHsl = kon_slh*ions_ca_sl*(Bmax_SLhighsl-buffers_SLHsl)-koff_slh*buffers_SLHsl;     // SLHsl     [mM/ms]
real J_CaB_junction = d_buffers_SLLj + d_buffers_SLHj;
real J_CaB_sl = d_buffers_SLLsl + d_buffers_SLHsl;

// Ion concentrations
// SR Ca Concentrations
real d_buffers_Csqn = kon_csqn*ions_ca_SR*(Bmax_Csqn-buffers_Csqn)-koff_csqn*buffers_Csqn;       // Csqn      [mM/ms]
real d_ions_ca_SR = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel) - d_buffers_Csqn;                     // Ca_sr      [mM/ms] %Ratio 3 leak current

// Sodium Concentrations
real I_Na_tot_junc = I_Na_junc+I_Nabk_junc+3*I_ncx_junc+3*I_nak_junc+ICaNa_junc;   // [uA/uF]
real I_Na_tot_sl = I_Na_sl+I_Nabk_sl+3*I_ncx_sl+3*I_nak_sl+ICaNa_sl;   // [uA/uF]

real d_ions_na_junc = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(ions_na_sl-ions_na_junc)-d_buffers_NaBj;
real d_ions_na_sl = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(ions_na_junc-ions_na_sl)+J_na_slmyo/Vsl*(ions_na_i-ions_na_sl)-d_buffers_NaBsl;
real d_ions_na_i = J_na_slmyo/Vmyo*(ions_na_sl-ions_na_i);             // [mM/msec]

// Potassium Concentration
real I_K_tot = I_to+IKr+I_ks+IK1-2*I_nak+I_CaK+I_Kb+Istim+IKatp;     // [uA/uF]
real d_ions_k_i = -I_K_tot*Cmem/(Vmyo*Frdy);           // [mM/msec]

// Cl concentrations
real I_Cl_tot = I_ClCa+I_Clbk;

real d_ions_cl_i = -I_Cl_tot*Cmem/(-1*Vmyo*Frdy);           // [mM/msec]

// Calcium Concentrations
real I_Ca_tot_junc = ICaL_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;  // [uA/uF]
real I_Ca_tot_sl = ICaL_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;   // [uA/uF]
real d_ions_ca_junc = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(ions_ca_sl-ions_ca_junc)-J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  // Ca_j
real d_ions_ca_sl = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(ions_ca_junc-ions_ca_sl) + J_ca_slmyo/Vsl*(ions_ca_i-ions_ca_sl)-J_CaB_sl;   // Ca_sl
real d_ions_ca_i = -J_serca*Vsr/Vmyo-J_CaB_cytosol+J_ca_slmyo/Vmyo*(ions_ca_sl-ions_ca_i);

// Membrane Potential
real I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          // [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk;                        // [uA/uF]
real I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
real I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
real d_membrane_v = -(I_tot);

// Derivative
rDY[0] = d_membrane_v;              // Euler
rDY[1] = d_ions_na_junc;            // Euler
rDY[2] = d_ions_na_sl;              // Euler
rDY[3] = d_ions_na_i;               // Euler
rDY[4] = d_ions_k_i;                // Euler
rDY[5] = d_ions_ca_junc;            // Euler
rDY[6] = d_ions_ca_sl;              // Euler
rDY[7] = d_ions_ca_i;               // Euler
rDY[8] = d_ions_cl_i;               // Euler
rDY[9] = d_ions_ca_SR;              // Euler
rDY[10] = d_buffers_NaBj;           // Euler  
rDY[11] = d_buffers_NaBsl;          // Euler
rDY[12] = d_buffers_TnClow;         // Euler
rDY[13] = d_buffers_TnCHc;          // Euler
rDY[14] = d_buffers_TnCHm;          // Euler
rDY[15] = d_buffers_CaM;            // Euler
rDY[16] = d_buffers_Myosin_ca;      // Euler
rDY[17] = d_buffers_Myosin_mg;      // Euler
rDY[18] = d_buffers_SRB;            // Euler
rDY[19] = d_buffers_SLLj;           // Euler
rDY[20] = d_buffers_SLLsl;          // Euler
rDY[21] = d_buffers_SLHj;           // Euler
rDY[22] = d_buffers_SLHsl;          // Euler
rDY[23] = d_buffers_Csqn;           // Euler
rDY[24] = dm;                       // RL
rDY[25] = dh;                       // RL
rDY[26] = dj;                       // RL
rDY[27] = dmL;                      // RL
rDY[28] = dhL;                      // RL
rDY[29] = dhLp;                     // RL
rDY[30] = d_ito_xtos;               // RL
rDY[31] = d_ito_ytos;               // RL
rDY[32] = d_ito_xtof;               // RL
rDY[33] = d_ito_ytof;               // RL
rDY[34] = d_ito_xtos_p;             // RL
rDY[35] = d_ito_ytos_p;             // RL
rDY[36] = d_ito_xtof_p;             // RL
rDY[37] = d_ito_ytof_p;             // RL
rDY[38] = d_iks_xs_junc;            // RL
rDY[39] = d_iks_xs_sl;              // RL
rDY[40] = dc0;                      // Euler
rDY[41] = dc1;                      // Euler
rDY[42] = dc2;                      // Euler
rDY[43] = delta_o;                  // Euler
rDY[44] = di;                       // Euler
rDY[45] = dd;                       // RL
rDY[46] = dff;                      // RL
rDY[47] = dfs;                      // RL
rDY[48] = dfcaf;                    // RL
rDY[49] = dfcas;                    // RL
rDY[50] = djca;                     // RL
rDY[51] = dnca;                     // Euler
rDY[52] = dnca_i;                   // Euler
rDY[53] = dffp;                     // RL
rDY[54] = dfcafp;                   // RL
rDY[55] = d_ical_pureCDI_junc;      // Euler
rDY[56] = d_ical_pureCDI_sl;        // Euler
rDY[57] = d_ryr_R;                  // Euler
rDY[58] = d_ryr_O;                  // Euler
rDY[59] = d_ryr_I;                  // Euler
rDY[60] = d_ryr_CaRI;               // Euler
rDY[61] = d_ryr_R_p;                // Euler
rDY[62] = d_ryr_O_p;                // Euler
rDY[63] = d_ryr_I_p;                // Euler
rDY[64] = d_ryr_CaRI_p;             // Euler
rDY[65] = d_jrel_icaldep_act;       // RL
rDY[66] = d_jrel_icaldep_f1;        // RL
rDY[67] = d_jrel_icaldep_f2;        // RL
rDY[68] = d_contraction_XS;         // Euler
rDY[69] = d_contraction_XW;         // Euler
rDY[70] = d_contraction_Ca_TRPN;    // Euler
rDY[71] = d_contraction_TmBlocked;  // Euler
rDY[72] = d_contraction_ZETAS;      // Euler
rDY[73] = d_contraction_ZETAW;      // Euler
rDY[74] = d_camk_trap;              // Euler
rDY[75] = d_camk_f_ICaL;            // RL
rDY[76] = d_camk_f_RyR;             // RL
rDY[77] = d_camk_f_PLB;             // RL
rDY[78] = d_casig_serca_trap;       // Euler
rDY[79] = dhp;                      // RL
rDY[80] = djp;                      // RL
rDY[81] = dm_P;                     // RL
rDY[82] = dh_P;                     // RL
rDY[83] = dj_P;                     // RL
rDY[84] = dhp_P;                    // RL
rDY[85] = djp_P;                    // RL
rDY[86] = dd_P;                     // RL
rDY[87] = dff_P;                    // RL
rDY[88] = dfs_P;                    // RL
rDY[89] = dfcaf_P;                  // RL
rDY[90] = dfcas_P;                  // RL
rDY[91] = dfBPf;                    // RL
rDY[92] = dfcaBPf;                  // RL