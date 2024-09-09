// Avoid NaN errors
//if (v >= 100.0) v = 100.0;
//if (v <= -100.0) v = -100.0;

//real aux, aux2, aux3;

// Changeable parameters
real nao = 140.0;
real cao = 1.8;
real ko = 5.0;
// Localization of ICaL and NCX: the fraction in junctional subspace
real ICaL_fractionSS = 0.8; 
real INaCa_fractionSS = 0.35;
// Additional parameters
real ca50 = 0.805;
real trpnmax = 0.07;

// INPUT CODE:
int mode = 0;           // 0 = "intact", 1 = "skinned"
real lambda = 1.0;
real lambda_rate = 0.0;

// EC parameters
real perm50 = 0.35;
real TRPN_n = 2;
real koff = 0.1;
real dr = 0.25;
real wfrac = 0.5;
real TOT_A = 25;
real ktm_unblock = 0.021; 
real beta_1 = -2.4;
real beta_0 = 2.3;
real gamma = 0.0085;
real gamma_wu = 0.615;
real phi = 2.23;

real nperm;
real Tref;
real nu;
real mu;
if (mode == 1) {
    nperm = 2.2;
    //ca50=2.5;
    Tref = 40.5;
    nu = 1;
    mu = 1;
}
else {
    nperm = 2.036; 
    //ca50=0.805;
    //ca50 = 0.5;
    //ca50 = 0.7;
    Tref = 120;
    nu = 7;
    mu = 3;
}

real k_ws = 0.004;
real k_uw = 0.026;

real lambda_min = 0.87;
real lambda_max = 1.2;

//real dydt[6] = {0,0,0,0,0,0};

k_ws = k_ws*mu;
k_uw = k_uw*nu;

real cdw = phi*k_uw*(1-dr)*(1-wfrac)/((1-dr)*wfrac);
real cds = phi*k_ws*(1-dr)*wfrac/dr;
real k_wu = k_uw*(1/wfrac-1)-k_ws;
real k_su = k_ws*(1/dr-1)*wfrac; 
real A = (0.25*TOT_A)/((1-dr)*wfrac+dr)*(dr/0.25);

//real lambda0 = (lambda_max < lambda) ? lambda_max : lambda;
real lambda0 = fminf(lambda_max, lambda);
//real lambda_aux = (lambda_min < lambda0) ? lambda_min : lambda0;
//aux = 1+beta_0*(lambda0+lambda_aux-(1+lambda_min));
//real Lfac = (aux > 0) ? aux : 0;
real Lfac = fmaxf(0,1+beta_0*(lambda0+fminf(lambda_min,lambda0)-(1+lambda_min)));

real XU = (1-TmBlocked)-XW-XS; // unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap
real xb_ws = k_ws*XW;
real xb_uw = k_uw*XU;
real xb_wu = k_wu*XW;
real xb_su = k_su*XS;

//aux = ((ZETAS>0)*ZETAS > (ZETAS<-1)*(-ZETAS-1)) ? (ZETAS>0)*ZETAS : (ZETAS<-1)*(-ZETAS-1); 
real gamma_rate=gamma*fmaxf((ZETAS>0)*ZETAS,(ZETAS<-1)*(-ZETAS-1));
real xb_su_gamma=gamma_rate*XS;
real gamma_rate_w=gamma_wu*abs(ZETAW); // weak xbs don't like being strained
real xb_wu_gamma=gamma_rate_w*XW;

//dydt[0] = xb_ws-xb_su-xb_su_gamma;
//dydt[1] = xb_uw-xb_wu-xb_ws-xb_wu_gamma;
real dXS = xb_ws-xb_su-xb_su_gamma;
real dXW = xb_uw-xb_wu-xb_ws-xb_wu_gamma;

//aux = (lambda-1 < 0.2) ? (lambda-1) : 0.2;
//ca50 = ca50+beta_1*aux;
ca50=ca50+beta_1*fminf(0.2,lambda-1);
//dydt[2] = koff*(pow(((cai*1000)/ca50),TRPN_n)*(1-Ca_TRPN)-Ca_TRPN); // untouched
real dCa_TRPN = koff*(pow(((cai*1000)/ca50),TRPN_n)*(1-Ca_TRPN)-Ca_TRPN);

real XSSS = dr*0.5;
real XWSS = (1-dr)*wfrac*0.5;
real ktm_block = ktm_unblock*(pow(perm50,nperm))*0.5/(0.5-XSSS-XWSS);

//aux = pow(Ca_TRPN,-(nperm/2));
//aux2 = pow(Ca_TRPN,(nperm/2));
//aux3 = (100 < aux) ? 100 : aux;
//dydt[3] = ktm_block*aux3*XU-ktm_unblock*aux2*TmBlocked;
//real dTmBlocked = ktm_block*aux3*XU-ktm_unblock*aux2*TmBlocked;
real dTmBlocked = ktm_block*fminf(100, pow(Ca_TRPN,-(nperm/2)))*XU - ktm_unblock*( pow(Ca_TRPN,(nperm/2)))*TmBlocked;

// velocity dependence -- assumes distortion resets on W->S
//dydt[4] = A*lambda_rate-cds*ZETAS;    // - gamma_rate * ZETAS;
//dydt[5] = A*lambda_rate-cdw*ZETAW;    // - gamma_rate_w * ZETAW;
real dZETAS = A*lambda_rate-cds*ZETAS;    // - gamma_rate * ZETAS;
real dZETAW = A*lambda_rate-cdw*ZETAW;    // - gamma_rate_w * ZETAW;

// Active Force
real Ta = Lfac*(Tref/dr)*((ZETAS+1)*XS+(ZETAW)*XW);

// physical constants
real R=8314.0;
real T=310.0;
real F=96485.0;

// cell geometry
real L=0.01;
real rad=0.0011;
real vcell=1000*3.14*rad*rad*L;
real Ageo=2*3.14*rad*rad+2*3.14*rad*L;
real Acap=2*Ageo;
real vmyo=0.68*vcell;
real vnsr=0.0552*vcell;
real vjsr=0.0048*vcell;
real vss=0.02*vcell;

real fkatp = 0.0;
real gkatp = 4.3195;

// CaMK constants
real KmCaMK=0.15;
real aCaMK=0.05*aCaMK_Multiplier;
real bCaMK=0.00068;
real CaMKo=0.05;
real KmCaM=0.0015;
// update CaMK
real CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
real CaMKa=CaMKb+CaMKt;
real dCaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt;      // Euler

// reversal potentials
real ENa=(R*T/F)*log(nao/nai);
real EK=(R*T/F)*log(ko/ki);
real PKNa=0.01833;
real EKs=(R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai));

// convenient shorthand calculations
real vffrt=v*F*F/(R*T);
real vfrt=v*F/(R*T);
real frt = F/(R*T);

real K_o_n = 5.0;
real A_atp = 2.0;
real K_atp = 0.25;
real akik = pow((ko / K_o_n), 0.24);
real bkik = (1.0 / (1.0 + pow((A_atp / K_atp), 2.0)));

real fINap=(1.0/(1.0+KmCaMK/CaMKa));
real fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
real fItop=(1.0/(1.0+KmCaMK/CaMKa));
real fICaLp=(1.0/(1.0+KmCaMK/CaMKa));

// INa formulations
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

// h phosphorylated
real hssp = 1 / pow((1 + exp( (v + 71.55 + 6)/7.43 )),2);
real dhp = (hssp - hp) / tauh;                  // Rush-Larsen
// j phosphorylated
real taujp = 1.46 * tauj;
real djp = (jss - jp) / taujp;                  // Rush-Larsen
real GNa = 11.7802;
real INa = INa_Multiplier*GNa*(v-ENa)*pow(m,3.0)*((1.0-fINap)*h*j+fINap*hp*jp);

// INaL
// calculate INaL
real mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
real tm = 0.1292 * exp(-pow(((v+45.79)/15.54),2)) + 0.06487 * exp(-pow(((v-4.823)/51.12),2)); 
real tmL=tm;
real dmL=(mLss-mL)/tmL;                                         // Rush-Larsen
real hLss=1.0/(1.0+exp((v+87.61)/7.488));
real thL=200.0;
real dhL=(hLss-hL)/thL;                                         // Rush-Larsen
real hLssp=1.0/(1.0+exp((v+93.81)/7.488));
real thLp=3.0*thL;
real dhLp=(hLssp-hLp)/thLp;                                     // Rush-Larsen
real GNaL=0.0279*INaL_Multiplier;
if (celltype==EPI) GNaL=GNaL*0.6;
real INaL = GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

// ITo
// calculate Ito
real ass=1.0/(1.0+exp((-(v-14.34))/14.82));
real ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
real da=(ass-a)/ta;                                                         // Rush-Larsen
real iss=1.0/(1.0+exp((v+43.94)/5.711));
real delta_epi = (celltype == EPI) ? 1.0-(0.95/(1.0+exp((v+70.0)/5.0))) : 1.0;
real tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
real tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
tiF=tiF*delta_epi;
tiS=tiS*delta_epi;
real AiF=1.0/(1.0+exp((v-213.6)/151.2));
real AiS=1.0-AiF;
real diF=(iss-iF)/tiF;                                                      // Rush-Larsen
real diS=(iss-iS)/tiS;                                                      // Rush-Larsen
real i=AiF*iF+AiS*iS;
real assp=1.0/(1.0+exp((-(v-24.34))/14.82));
real dap=(assp-ap)/ta;                                                     // Rush-Larsen
real dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
real dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
real tiFp=dti_develop*dti_recover*tiF;
real tiSp=dti_develop*dti_recover*tiS;
real diFp=(iss-iFp)/tiFp;                                                   // Rush-Larsen
real diSp=(iss-iSp)/tiSp;                                                   // Rush-Larsen
real ip=AiF*iFp+AiS*iSp;
real Gto=0.16*Ito_Multiplier;
Gto = (celltype == EPI || celltype == MID) ? Gto*2.0 : Gto;
real Ito = Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

// ICaL
// a variant updated by jakub, using a changed activation curve
// it computes both ICaL in subspace and myoplasm (_i)

// calculate ICaL, ICaNa, ICaK
real dss=1.0763*exp(-1.0070*exp(-0.0829*(v)));  // magyar
if(v >31.4978) dss = 1; // activation cannot be greater than 1
real td= 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
real dd=(dss-d)/td;                                                                     // Rush-Larsen
real fss=1.0/(1.0+exp((v+19.58)/3.696));
real tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
real tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
real Aff=0.6;
real Afs=1.0-Aff;
real dff=(fss-ff)/tff;                                                                  // Rush-Larsen
real dfs=(fss-fs)/tfs;                                                                  // Rush-Larsen
real f=Aff*ff+Afs*fs;
real fcass=fss;
real tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
real tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));

real Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
real Afcas=1.0-Afcaf;
real dfcaf=(fcass-fcaf)/tfcaf;                                                          // Rush-Larsen
real dfcas=(fcass-fcas)/tfcas;                                                          // Rush-Larsen
real fca=Afcaf*fcaf+Afcas*fcas;

real tjca = 75;
//real tjca = 72.5;
real jcass = 1.0/(1.0+exp((v+18.08)/(2.7916)));   
real djca=(jcass-jca)/tjca;                                                                  // Rush-Larsen
real tffp=2.5*tff;
real dffp=(fss-ffp)/tffp;                                                                    // Rush-Larsen
real fp=Aff*ffp+Afs*fs;
real tfcafp=2.5*tfcaf;
real dfcafp=(fcass-fcafp)/tfcafp;                                                           // Rush-Larsen
real fcap=Afcaf*fcafp+Afcas*fcas;

// SS nca
real Kmn=0.002;
real k2n=500.0;
real km2n=jca*1;
real anca=1.0/(k2n/km2n+pow((1.0+Kmn/cass),4.0));
real dnca=anca*k2n-nca*km2n;                                                        // Euler

// myoplasmic nca
real anca_i = 1.0/(k2n/km2n+pow((1.0+Kmn/cai),4.0));
real dnca_i = anca_i*k2n-nca_i*km2n;                                                     // Euler

// SS driving force
//real clo = 150;     // Extracellular Cl  [mM]
//real cli = 24.0;    // Intracellular Cl  [mM]
real Io = 0.5*(nao + ko + clo + 4*cao)/1000;         // ionic strength outside. /1000 is for things being in micromolar
real Ii = 0.5*(nass + kss + cli + 4*cass)/1000;     // ionic strength outside. /1000 is for things being in micromolar
//real Ii = 0.5*(nass + kss + clss + 4*cass)/1000;     // (dynCl) ionic strength outside. /1000 is for things being in micromolar
// The ionic strength is too high for basic DebHuc. We'll use Davies
real dielConstant = 74;     // water at 37 degrees.
real temp = 310;            // body temp in kelvins.
real constA = 1.82*pow(10,6)*pow((dielConstant*temp),(-1.5));

real gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
real gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
real gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
real gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
real gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
real gamma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

real PhiCaL_ss =  4.0*vffrt*(gamma_cai*cass*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
real PhiCaNa_ss =  1.0*vffrt*(gamma_nai*nass*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
real PhiCaK_ss =  1.0*vffrt*(gamma_ki*kss*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);

// Myo driving force
Io = 0.5*(nao + ko + clo + 4*cao)/1000; // ionic strength outside. /1000 is for things being in micromolar
Ii = 0.5*(nai + ki + cli + 4*cai)/1000; // ionic strength outside. /1000 is for things being in micromolar
// The ionic strength is too high for basic DebHuc. We'll use Davies

gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

real gammaCaoMyo = gamma_cao;
real gammaCaiMyo = gamma_cai;

real PhiCaL_i =  4.0*vffrt*(gamma_cai*cai*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
real PhiCaNa_i =  1.0*vffrt*(gamma_nai*nai*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
real PhiCaK_i =  1.0*vffrt*(gamma_ki*ki*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);
// The rest
real PCa=8.3757e-05 * ICaL_Multiplier;
if (celltype==EPI)
    PCa=PCa*1.2;
else if (celltype==MID)
    //PCa=PCa*2;
    PCa = PCa*1.8;

real PCap=1.1*PCa;
real PCaNa=0.00125*PCa;
real PCaK=3.574e-4*PCa;
real PCaNap=0.00125*PCap;
real PCaKp=3.574e-4*PCap;

real ICaL_ss=(1.0-fICaLp)*PCa*PhiCaL_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
real ICaNa_ss=(1.0-fICaLp)*PCaNa*PhiCaNa_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
real ICaK_ss=(1.0-fICaLp)*PCaK*PhiCaK_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK_ss*d*(fp*(1.0-nca)+jca*fcap*nca);

real ICaL_i=(1.0-fICaLp)*PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
real ICaNa_i=(1.0-fICaLp)*PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
real ICaK_i=(1.0-fICaLp)*PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);

// And we weight ICaL (in ss) and ICaL_i
ICaL_i = ICaL_i * (1-ICaL_fractionSS);
ICaNa_i = ICaNa_i * (1-ICaL_fractionSS);
ICaK_i = ICaK_i * (1-ICaL_fractionSS);
ICaL_ss = ICaL_ss * ICaL_fractionSS;
ICaNa_ss = ICaNa_ss * ICaL_fractionSS;
ICaK_ss = ICaK_ss * ICaL_fractionSS;

real ICaL = ICaL_ss + ICaL_i;
real ICaNa = ICaNa_ss + ICaNa_i;
real ICaK = ICaK_ss + ICaK_i;
real ICaL_tot = ICaL + ICaNa + ICaK;

// Ikr
// Variant based on Lu-Vandenberg
// Extracting state vector
// IMPORTANT: For reason of backward compatibility of naming of an older version of a MM IKr
//            c3 in code is c0 in article diagram, c2 is c1, c1 is c2
real c0 = ikr_c0;
real c1 = ikr_c1;
real c2 = ikr_c2;
//real c0 = ikr_c2;
//real c1 = ikr_c1;
//real c2 = ikr_c0;
real o = ikr_o;
real I = ikr_i;
real b = 0; // no channels blocked in via the mechanism of specific MM states

// transition rates
// from c0 to c1 in l-v model,
real alpha = 0.1161 * exp(0.2990 * vfrt);
// from c1 to c0 in l-v/
real beta =  0.2442 * exp(-1.604 * vfrt);

// from c1 to c2 in l-v/
real alpha1 = 1.25 * 0.1235;
// from c2 to c1 in l-v/
real beta1 =  0.1911;

// from c2 to o/           c1 to o
real alpha2 =0.0578 * exp(0.9710 * vfrt);
// from o to c2/
real beta2 = 0.349e-3* exp(-1.062 * vfrt);

// from o to i
real alphai = 0.2533 * exp(0.5953 * vfrt);
// from i to o
real betai = 1.25* 0.0522 * exp(-0.8209 * vfrt);

// from c2 to i (from c1 in orig)
real alphac2ToI = 0.52e-4 * exp(1.525 * vfrt);
// from i to c2
real betaItoC2 = 0.85e-8 * exp(-1.842 * vfrt);
betaItoC2 = (beta2 * betai * alphac2ToI)/(alpha2 * alphai);
// transitions themselves
// for reason of backward compatibility of naming of an older version of a
// MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

real dc0 = c1 * beta - c0 * alpha;                   // Euler
real dc1 = c0 * alpha + c2*beta1 - c1*(beta+alpha1); // Euler
real dc2 = c1 * alpha1 + o*beta2 + I*betaItoC2 - c2 * (beta1 + alpha2 + alphac2ToI); // Euler
real delta_o = c2 * alpha2 + I*betai - o*(beta2+alphai);    // Euler
real di = c2*alphac2ToI + o*alphai - I*(betaItoC2 + betai); // Euler

real GKr = 0.0321 * sqrt(ko/5) * IKr_Multiplier; // 1st element compensates for change to ko (sqrt(5/5.4)* 0.0362)
if (celltype==EPI)
    GKr=GKr*1.3;
else if (celltype==MID)
    GKr=GKr*0.8;

real IKr = GKr * o  * (v-EK);

// calculate IKs
real xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
real txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
real dxs1=(xs1ss-xs1)/txs1;                              // Rush-Larsen
real xs2ss=xs1ss;
real txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
real dxs2=(xs2ss-xs2)/txs2;                               // Rush-Larsen
real KsCa=1.0+0.6/(1.0+pow((3.8e-5/cai),1.4));
real GKs= 0.0011*IKs_Multiplier;
if (celltype==EPI)
    GKs=GKs*1.4;
real IKs = GKs*KsCa*xs1*xs2*(v-EKs);

// IK1
real aK1 = 4.094/(1+exp(0.1217*(v-EK-49.934)));
real bK1 = (15.72*exp(0.0674*(v-EK-3.257))+exp(0.0618*(v-EK-594.31)))/(1+exp(-0.1629*(v-EK+14.207)));
real K1ss = aK1/(aK1+bK1);
real GK1=IK1_Multiplier  * 0.6992; // 0.7266; // * sqrt(5/5.4))
if (celltype==EPI)
    GK1=GK1*1.2;
else if (celltype==MID)
    GK1=GK1*1.3;
real IK1=GK1*sqrt(ko/5)*K1ss*(v-EK);

// IKCa
real fIKCass = 0.8;
real kdikca = 6.05e-4;
real ikcan = 3.5;
real GKCa = 0.003 * IKCa_Multiplier;
real IKCa_ss = GKCa * fIKCass * pow(cass,ikcan) / (pow(cass,ikcan) + pow(kdikca,ikcan)) * (v-EK);
real IKCa_i = GKCa * (1.0-fIKCass) * pow(cai,ikcan) / (pow(cai,ikcan) + pow(kdikca,ikcan)) * (v-EK);
real IKCa = IKCa_ss + IKCa_i;

// INaCa
real zca = 2.0;
real kna1=15.0;
real kna2=5.0;
real kna3=88.12;
real kasymm=12.5;
real wna=6.0e4;
real wca=6.0e4;
real wnaca=5.0e3;
real kcaon=1.5e6;
real kcaoff=5.0e3;
real qna=0.5224;
real qca=0.1670;
real hca=exp((qca*v*F)/(R*T));
real hna=exp((qna*v*F)/(R*T));

real h1=1+nai/kna3*(1+hna);
real h2=(nai*hna)/(kna3*h1);
real h3=1.0/h1;
real h4=1.0+nai/kna1*(1+nai/kna2);
real h5=nai*nai/(h4*kna1*kna2);
real h6=1.0/h4;
real h7=1.0+nao/kna3*(1.0+1.0/hna);
real h8=nao/(kna3*hna*h7);
real h9=1.0/h7;
real h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
real h11=nao*nao/(h10*kna1*kna2);
real h12=1.0/h10;
real k1=h12*cao*kcaon;
real k2=kcaoff;
real k3p=h9*wca;
real k3pp=h8*wnaca;
real k3=k3p+k3pp;
real k4p=h3*wca/hca;
real k4pp=h2*wnaca;
real k4=k4p+k4pp;
real k5=kcaoff;
real k6=h6*cai*kcaon;
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
real KmCaAct=150.0e-6;
real allo=1.0/(1.0+pow((KmCaAct/cai),2.0));
real zna=1.0;
real JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
real JncxCa=E2*k2-E1*k1;
real Gncx= 0.0034 * INaCa_Multiplier;
if (celltype==EPI)
    Gncx=Gncx*1.1;
else if (celltype==MID)
    Gncx=Gncx*1.4;
real INaCa_i = (1-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa);

// calculate INaCa_ss
h1=1+nass/kna3*(1+hna);
h2=(nass*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+nass/kna1*(1+nass/kna2);
h5=nass*nass/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+nao/kna3*(1.0+1.0/hna);
h8=nao/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
h11=nao*nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*cao*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*cass*kcaon;
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
KmCaAct=150.0e-6 ;
allo=1.0/(1.0+pow((KmCaAct/cass),2.0));
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
real INaCa_ss = INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa);

// calculate INaK
real k1p=949.5;
real k1m=182.4;
real k2p=687.2;
real k2m=39.4;
k3p=1899.0;
real k3m=79300.0;
k4p=639.0;
real k4m=40.0;
real Knai0=9.073;
real Knao0=27.78;
real delta=-0.1550;
real Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
real Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
real Kki=0.5;
real Kko=0.3582;
real MgADP=0.05;
real MgATP=9.8;
real Kmgatp=1.698e-7;
real H=1.0e-7;
real eP=4.2;
real Khp=1.698e-7;
real Knap=224.0;
real Kxkur=292.0;
real P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
real a1=(k1p*pow((nai/Knai),3.0))/(pow((1.0+nai/Knai),3.0)+pow((1.0+ki/Kki),2.0)-1.0);
real b1=k1m*MgADP;
real a2=k2p;
real b2=(k2m*pow((nao/Knao),3.0))/(pow((1.0+nao/Knao),3.0)+pow((1.0+ko/Kko),2.0)-1.0);
real a3=(k3p*pow((ko/Kko),2.0))/(pow((1.0+nao/Knao),3.0)+pow((1.0+ko/Kko),2.0)-1.0);
real b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
real a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
real b4=(k4m*pow((ki/Kki),2.0))/(pow((1.0+nai/Knai),3.0)+pow((1.0+ki/Kki),2.0)-1.0);
x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
real zk=1.0;
real JnakNa=3.0*(E1*a3-E2*b3);
real JnakK=2.0*(E4*b1-E3*a1);
real Pnak= 15.4509;
if (celltype==EPI)
    Pnak=Pnak*0.9;
else if (celltype==MID)
    Pnak=Pnak*0.7;

real INaK = Pnak*(zna*JnakNa+zk*JnakK)*INaK_Multiplier;

// Minor/background currents
// calculate IKb
real xkb=1.0/(1.0+exp(-(v-10.8968)/(23.9871)));
real GKb=0.0189*IKb_Multiplier;

if (IKCa_Multiplier > 0.0)
    GKb = GKb*0.9;
if (celltype==EPI)
    GKb=GKb*0.6;
real IKb = GKb*xkb*(v-EK);

// calculate INab
real PNab=1.9239e-09*INab_Multiplier;
real INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

// calculate ICab
real PCab=5.9194e-08*ICab_Multiplier; 
real ICab=PCab*4.0*vffrt*(gammaCaiMyo*cai*exp(2.0*vfrt)-gammaCaoMyo*cao)/(exp(2.0*vfrt)-1.0);

// calculate IpCa
real GpCa=5e-04*IpCa_Multiplier;
real IpCa=GpCa*cai/(0.0005+cai);

// Chloride
// I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current

real ECl = (R*T/F)*log(cli/clo);            // [mV]
//real EClss = (R*T/F)*log(clss/clo);         // [mV]   dynCl

real Fjunc = 1;   
real Fsl = 1-Fjunc; // fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace

real GClCa = ICaCl_Multiplier * 0.2843;   // [mS/uF]
real GClB = IClb_Multiplier * 1.98e-3;    // [mS/uF] 
real KdClCa = 0.1;                        // [mM]

//real I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/cass)*(v-EClss);     // dynCl
//real I_ClCa_sl = Fsl*GClCa/(1+KdClCa/cai)*(v-ECl);            // dynCl
real I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/cass)*(v-ECl);
real I_ClCa_sl = Fsl*GClCa/(1+KdClCa/cai)*(v-ECl);

real I_ClCa = I_ClCa_junc+I_ClCa_sl;
real I_Clbk = GClB*(v-ECl);

// Calcium handling
// calculate ryanodione receptor calcium induced calcium release from the jsr
real fJrelp=(1.0/(1.0+KmCaMK/CaMKa));

// Jrel
real jsrMidpoint = 1.7;

real bt=4.75;
real a_rel=0.5*bt;
real Jrel_inf=a_rel*(-ICaL)/(1.0+pow((jsrMidpoint/cajsr),8.0));
if (celltype==MID)
    Jrel_inf=Jrel_inf*1.7;
real tau_rel=bt/(1.0+0.0123/cajsr);

if (tau_rel<0.001)
    tau_rel=0.001;

real dJrelnp=(Jrel_inf-Jrel_np)/tau_rel;                     // Rush-Larsen

real btp=1.25*bt;
real a_relp=0.5*btp;
real Jrel_infp=a_relp*(-ICaL)/(1.0+pow((jsrMidpoint/cajsr),8.0));
if (celltype==MID)
    Jrel_infp=Jrel_infp*1.7;
real tau_relp=btp/(1.0+0.0123/cajsr);

if (tau_relp<0.001)
    tau_relp=0.001;

tau_relp = tau_relp*taurelp_Multiplier;
real dJrelp=(Jrel_infp-Jrel_p)/tau_relp;                     // Rush-Larsen
real Jrel = Jrel_Multiplier * 1.5378 * ((1.0-fJrelp)*Jrel_np+fJrelp*Jrel_p);

real fJupp=(1.0/(1.0+KmCaMK/CaMKa));

// Jup
// calculate serca pump, ca uptake flux
// camkFactor = 2.4;
// gjup = 0.00696;
// Jupnp=Jup_Multiplier * gjup*cai/(cai+0.001);
// Jupp=Jup_Multiplier * camkFactor*gjup*cai/(cai + 8.2500e-04);
// if celltype==1
//     Jupnp=Jupnp*1.3;
//     Jupp=Jupp*1.3;
// end
// 
// 
// Jleak=Jup_Multiplier * 0.00629 * cansr/15.0;
// Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

// calculate serca pump, ca uptake flux
real Jupnp = Jup_Multiplier * 0.005425*cai/(cai+0.00092);
real Jupp = Jup_Multiplier * 2.75*0.005425*cai/(cai+0.00092-0.00017);
if (celltype==EPI) {
    Jupnp=Jupnp*1.3;
    Jupp=Jupp*1.3;
}
real Jleak=Jup_Multiplier* 0.0048825*cansr/15.0;
real Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

//calculate tranlocation flux
real Jtr=(cansr-cajsr)/60;

// I_katp current (fkatp current)
//real I_katp = (fkatp * gkatp * akik * bkik * (v - EK));

// stimulus current
real Istim = calc_I_stim;

//update the membrane voltage
//real dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+Istim);    // Euler
real dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+IKCa+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+Istim);    // Euler
//real dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+I_katp+Istim);    // Euler

// calculate diffusion fluxes
real JdiffNa=(nass-nai)/2.0;
real JdiffK=(kss-ki)/2.0;
real Jdiff=(cass-cai)/0.2;
//real JdiffCl=(clss-cli)/2.0;  // dynCl

// calcium buffer constants 
real cmdnmax= 0.05; 
if (celltype==EPI)
    cmdnmax=cmdnmax*1.3;
real kmcmdn=0.00238; 
//real trpnmax=0.07;
real kmtrpn=0.0005;
real BSRmax=0.047;
real KmBSR = 0.00087;
real BSLmax=1.124;
real KmBSL = 0.0087;
real csqnmax=10.0;
real kmcsqn=0.8;

// update intracellular concentrations, using buffers for cai, cass, cajsr
real dnai=-(ICaNa_i+INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;     // Euler
real dnass=-(ICaNa_ss+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;                                   // Euler

real dki=-(ICaK_i+Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;        // Euler
real dkss=-(ICaK_ss)*Acap/(F*vss)-JdiffK;                                                   // Euler

//real Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow((kmcmdn+cai),2.0)+trpnmax*kmtrpn/pow((kmtrpn+cai),2.0)); // dynCl
real Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow((kmcmdn+cai),2.0));
//real dcai=Bcai*(-(ICaL_i + IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo); // dynCl
real dcai=Bcai*(-(ICaL_i + IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo-dCa_TRPN*trpnmax);

real Bcass=1.0/(1.0+BSRmax*KmBSR/pow((KmBSR+cass),2.0)+BSLmax*KmBSL/pow((KmBSL+cass),2.0));
real dcass=Bcass*(-(ICaL_ss-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);

real dcansr=Jup-Jtr*vjsr/vnsr;

real Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow((kmcsqn+cajsr),2.0));
real dcajsr=Bcajsr*(Jtr-Jrel);

//real dcli  = - (I_Clbk + I_ClCa_sl)*Acap/(-1*F*vmyo)+JdiffCl*vss/vmyo;        // dynCl
//real dclss = - I_ClCa_junc*Acap/(-1*F*vss)-JdiffCl;                           // dynCl

// Right-hand side
rDY_[0]  = dv;
rDY_[1]  = dnai;
rDY_[2]  = dnass;
rDY_[3]  = dki;
rDY_[4]  = dkss;
rDY_[5]  = dcai;
rDY_[6]  = dcass;
rDY_[7]  = dcansr;
rDY_[8]  = dcajsr;
rDY_[9]  = dm;
rDY_[10] = dhp;
rDY_[11] = dh;
rDY_[12] = dj;
rDY_[13] = djp;
rDY_[14] = dmL;
rDY_[15] = dhL;
rDY_[16] = dhLp;
rDY_[17] = da;
rDY_[18] = diF;
rDY_[19] = diS;
rDY_[20] = dap;
rDY_[21] = diFp;
rDY_[22] = diSp;
rDY_[23] = dd;
rDY_[24] = dff;
rDY_[25] = dfs;
rDY_[26] = dfcaf;
rDY_[27] = dfcas;
rDY_[28] = djca;
rDY_[29] = dnca;
rDY_[30] = dnca_i;
rDY_[31] = dffp;
rDY_[32] = dfcafp;
rDY_[33] = dxs1;
rDY_[34] = dxs2;
rDY_[35] = dJrelnp;
rDY_[36] = dCaMKt;
rDY_[37] = dc0;
rDY_[38] = dc1;
rDY_[39] = dc2;
rDY_[40] = delta_o;
rDY_[41] = di;
rDY_[42] = dJrelp;
// -----------------------
// Land-Niederer
rDY_[43] = dXS;
rDY_[44] = dXW;
rDY_[45] = dCa_TRPN;
rDY_[46] = dTmBlocked;
rDY_[47] = dZETAS;
rDY_[48] = dZETAW;
rDY_[49] = Ta;
