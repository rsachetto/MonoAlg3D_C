// Avoid NaN errors
if (v >= 100.0) v = 100.0;
if (v <= -100.0) v = -100.0;

real INaCa_fractionSS = 0.2;

// Constants
real pi = 3.14;
real cao = 1.8;  //[Ca]o mM
real nao = 140.0;  //[Na]o mM
real ko  = 5.4;  //[K]o  mM
// Physical Constants:
real R = 8314.0;   // J/kmol/K
real T = 310.0;    // K
real F = 96485.0;  // C/mol
real vffrt = v*F*F/(R*T);
real vfrt  = v*F/(R*T);
// Valence of ions
real zna = 1;
real zk  = 1;
real zca = 2;
// Cell Geometry from ORd
// Cell geometry was approxymate by a cylinder of length L and radius r
real L = 0.0164;                         // cm
real rad = 0.00175;                      // cm
real vcell = 1000*pi*rad*rad*L;            // uL
// Geometric Area
real Ageo = 2*pi*rad*rad + 2*pi*rad*L;     // cm^2
// Capacitive Area
real Acap = 2*Ageo;                      // cm^2
// Compartment Volumes (6)
real vmyo = 0.60*vcell;                  // uL
real vnsr = 0.04*vcell;                  // uL
real vjsr = 0.002*vcell;                 // uL
real vcsr = 0.008*vcell;                   // uL
real vss  = 0.02*vcell;                  // uL
real vsl  = 0.15*vcell;                     // uL
// Reversal Potentials
real ENa  = (R*T/F)*log(nao/nasl);
real EK   = (R*T/F)*log(ko/ksl);
real ECa= (R*T/(zca*F))*log(cao/casl);
real PKNa = 0.01833;
real EKs  = (R*T/F)*log((ko+PKNa*nao)/(ksl+PKNa*nasl));
// Calcium_Fluxes_rate_constants from PRd
real tautr1 = 120;
real tautr2 = 120;	
real gaptau = 12;
real sstau = 0.2;
// CaMK Constants
real KmCaMK = 0.15;  
real aCaMK  = 0.05;  
real bCaMK  = 0.00068;
real CaMKo  = 0.05;  
real KmCaM  = 0.0015;
// update CaMK -> X(41)
real CaMKb = CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
real CaMKa = CaMKb+CaMKt;
real dCaMKt = aCaMK*CaMKb*(CaMKb+CaMKt) - bCaMK*CaMKt;

// ------------------------------------------------------------------------
//                      MEMBRANE IONIC CURRENTS                        
// ------------------------------------------------------------------------
// INa current from ORd
real mss=1.0/ (1.0+exp((-(v+48.4264))/7.5653)); 
real tm=1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
real hss=1.0/(1+exp((v+78.5)/6.22)); 
real thf=1.0/(3.6860e-06*exp(-(v+3.8875)/7.8579)+16*exp((v-0.4963)/9.1843));
real ths=1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
real Ahf=0.99;
real Ahs=1.0-Ahf;
real h=Ahf*hf+Ahs*hs;
real jss=hss;
real tj=4.8590+1.0/(0.8628*exp(-(v+116.7258)/7.6005)+1.1096*exp((v+6.2719)/9.0358)); 
real hssp=1.0/(1+exp((v+84.7)/6.22));
real thsp=3.0*ths;
real hp=Ahf*hf+Ahs*hsp;
real tjp=1.46*tj;
real dm=(mss-m)/tm;
real dhf=(hss-hf)/thf;
real dhs=(hss-hs)/ths;
real dj=(jss-j)/tj;
real dhsp=(hssp-hsp)/thsp;
real djp=(jss-jp)/tjp;
real GNa=39.4572*GNa_Multiplier;
real fINap=(1.0/(1.0+KmCaMK/CaMKa));
real INa=GNa*(v-ENa)*pow(m,3.0)*((1.0-fINap)*h*j+fINap*hp*jp);
INa *= INa_Multiplier;

// INaL current [from ORd]
real mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
real tmL=tm;
real hLss=1.0/(1.0+exp((v+87.61)/7.488));
real thL=200.0;
real hLssp=1.0/(1.0+exp((v+93.81)/7.488));
real thLp=3.0*thL;
real dmL=(mLss-mL)/tmL;
real dhL=(hLss-hL)/thL;
real dhLp=(hLssp-hLp)/thLp;
real GNaL=0.0189*GNaL_Multiplier; 
real fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
real INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
INaL *= INaL_Multiplier;

// ICaL [from ORd]
real vshift=15.19;
real dss=1.0/(1.0+exp((-(v+3.940+3.3))/4.230));
real td=0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
real fss=1.0/(1.0+exp((v+19.58+3.3)/3.696));
real tff=7.0+1.0/(0.0045*exp(-((v+vshift)+20.0)/10.0)+0.0045*exp(((v+vshift)+20.0)/10.0));
real tfs=1000.0+1.0/(0.000035*exp(-((v+vshift)+5.0)/4.0)+0.000035*exp(((v+vshift)+5.0)/6.0)); 
real Aff=0.6;
real Afs=1.0-Aff;
real f=Aff*ff+Afs*fs;
real fcass=fss;
real tfcaf=0.72*(7.0+1.0/(0.04*exp(-((v+vshift)-4.0)/7.0)+0.04*exp(((v+vshift)-4.0)/7.0)));
real tfcas=0.49*(100.0 + 1.0/(0.00012*exp(-(v+vshift)/3.0)+0.00012*exp((v+vshift)/7.0)));
real Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
real Afcas=1.0-Afcaf;
real fca=Afcaf*fcaf+Afcas*fcas;
real tjca=75.0;
real ktaup=2.5;
real tffp=ktaup*tff;
real fp=Aff*ffp+Afs*fs;
real tfcafp=ktaup*tfcaf;
real fcap=Afcaf*fcafp+Afcas*fcas;
real Kmn=0.002;
real k2n=1000.0;
real km2n=jca*1.0;
real anca=1.0/(k2n/km2n+pow((1.0+Kmn/cass),4.0));
real dnca=(anca*k2n-nca*km2n);
real PhiCaL=4.0*vffrt*(cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
real PhiCaNa=1.0*vffrt*(0.75*nass*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
real PhiCaK=1.0*vffrt*(0.75*kss*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
real PCa=0.0001*0.776772893145948;
real PCap=1.1*PCa;
real PCaNa=0.00125*PCa;
real PCaK=3.574e-4*PCa;
real PCaNap=0.00125*PCap;
real PCaKp=3.574e-4*PCap;
real dd=(dss-d)/td;
real dff=(fss-ff)/tff;
real dfs=(fss-fs)/tfs;
real dfcaf=(fcass-fcaf)/tfcaf;
real dfcas=(fcass-fcas)/tfcas;
real djca=(fcass-jca)/tjca;
real dffp=(fss-ffp)/tffp;
real dfcafp=(fcass-fcafp)/tfcafp;
real fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
real ICaL=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
real ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
real ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaL *= ICaL_Multiplier;
ICaNa *= ICaNa_Multiplier;
ICaK *= ICaK_Multiplier;

// T_type_calcium_current_ICaT [from PRd]
real gcat= 0.0754*GCaT_Multiplier;      
real bss= 1/(1+ exp (-(v+30)/7));
real gss= 1/(1+exp((v+61)/5));
real taub= 1/(1.068*exp((v+16.3)/30)+1.068*exp(-(v+16.3)/30));
real taug= 1/(0.015*exp(-(v+71.7)/83.3)+0.015*exp((v+71.7)/15.4));
real db= (bss-b)/taub;
real dg= (gss-g)/taug;
real ICaT= gcat*b*g*(v-ECa);
ICaT *= ICaT_Multiplier;

// Ito current [from Han's data]
real gtos= 0.192*Gto_Multiplier;          
real ass= 1.0/(1.0+exp((20.0-v)/13.0));
real iss= 1.0/(1.0+exp((v+27.0)/13.0));
real atau= 1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+ 3.5/(1.0+exp((v+100.0)/29.3814))); 
real tiS=43+1./(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
real tiF=6.162+1./(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v-8.0)/8.59));
real da= (ass-a)/atau;
real di1= (iss-i1)/tiS;
real di2= (iss-i2)/tiF;
real Ito= gtos*a*i1*i2*(v-EK);
Ito *= Ito_Multiplier;

// Sustained_outward_current [from Han's data]
real g_sus=0.0301*Gsus_Multiplier;
real asus = 1.0/(1.0+exp(-(v-12)/16.0));
real Isus= g_sus*asus*(v-EK);
Isus *= Isus_Multiplier;

// IKr current [modified from ORd]
real xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
real txrf=12.98+1.0/(0.3652*exp((v-31.66+17.6)/3.869)+4.123e-5*exp((-(v-47.78+17.6))/20.38));
real txrs=1.865+1.0/(0.06629*exp((v-34.70+17.2)/7.355)+1.128e-5*exp((-(v-29.74+17.2))/25.94));
real Axrf=1.0/(1.0+exp((v+54.81)/38.21));
real Axrs=1.0-Axrf;
real dxrf=(xrss-xrf)/txrf;
real dxrs=(xrss-xrs)/txrs;
real xr=Axrf*xrf+Axrs*xrs;
real rkr=1/(1+exp((v+55)/(0.32*75)))*1/(1+exp((v-10)/(0.32*30)));
real GKr=0.0342*Gkr_Multiplier;
real IKr=GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);
IKr *= IKr_Multiplier;

// IKs current [from ORd]
real xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
real txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
real dxs1=(xs1ss-xs1)/txs1;
real xs2ss=xs1ss;
real txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
real dxs2=(xs2ss-xs2)/txs2;
real KsCa=1.0+0.6/(1.0+pow(3.8e-5/casl,1.4));
real GKs=0.0029*Gks_Multiplier;                       
real IKs=GKs*KsCa*xs1*xs2*(v-EKs);
IKs *= IKs_Multiplier;

// Hyper-polarization Activated Current (If) [From PRd]
real yss= 1/(1+exp((v+87)/9.5));
real tau_y= 2000/(exp(-(v+132)/10) + exp((v+57)/60));
real dy=(yss-y)/tau_y;
real g_fna = 0.0116;
g_fna *= GfNa_Multiplier;
real g_fk = 0.0232; 
g_fk *= GfK_Multiplier;
real IfNa  = g_fna*y*y*(v-ENa);   
real IfK = g_fk*y*y*(v-EK);    
real If = IfNa + IfK;
If *= If_Multiplier;

// IK1 current [from ORd's formulation with Han's data]
real gk1=2.3238*sqrt(ko/5.4)*0.0455;      
gk1 *= GK1_Multiplier;         
real xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
real xk1tau=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
real dxk1=(xk1ss-xk1)/xk1tau;
real rk1=1/(1.0+exp((v+116-5.5*ko)/11)); 
real IK1=gk1*rk1*xk1*(v-EK);
IK1 *= IK1_Multiplier;

// INaCa current [from ORd]
real kna1=15.0;      
real kna2=5.0;       
real kna3=88.12;     
real kasymm=12.5;
real wna=6.0e4;      
real wca=6.0e4;      
real wnaca=5.0e3;    
real KmCaAct=150.0e-6;
real kcaon=1.5e6;    
real kcaoff=5.0e3;   
real qna=0.5224;     
real qca=0.1670;
real Gncx=9.5709e-04;
Gncx *= GNCX_Multiplier;

// INaCa_i current
real hca=exp((qca*v*F)/(R*T));       
real hna=exp((qna*v*F)/(R*T));
real h1=1+nasl/kna3*(1+hna);          
real h2=(nasl*hna)/(kna3*h1);
real h3=1.0/h1;                      
real h4=1.0+nasl/kna1*(1+nasl/kna2);
real h5=nasl*nasl/(h4*kna1*kna2);      
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
real k6=h6*casl*kcaon;  
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

real allo=1.0/(1.0+pow((KmCaAct/casl),2.0));   
real JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;   
real JncxCa=E2*k2-E1*k1;
real INaCa_i=(1-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa);

// INaCa_ss current
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

allo=1.0/(1.0+pow((KmCaAct/cass),2.0));
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
real INaCa_ss=INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa);
real INaCa = INaCa_i + INaCa_ss;
INaCa *= INaCa_Multiplier;

// INaK current [from ORd]
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
real delta2=-0.1550;  
real Knai=Knai0*exp((delta2*v*F)/(3.0*R*T));
real Knao=Knao0*exp(((1.0-delta2)*v*F)/(3.0*R*T));
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
real P=eP/(1.0+H/Khp+nasl/Knap+ksl/Kxkur);

real a1=(k1p*pow((nasl/Knai),3.0))/(pow((1.0+nasl/Knai),3.0)+pow((1.0+ksl/Kki),2.0)-1.0);
real b1=k1m*MgADP;
real a2=k2p;
real b2=(k2m*pow((nao/Knao),3.0))/(pow((1.0+nao/Knao),3.0)+pow((1.0+ko/Kko),2.0)-1.0);
real a3=(k3p*pow((ko/Kko),2.0))/(pow((1.0+nao/Knao),3.0)+pow((1.0+ko/Kko),2.0)-1.0);
real b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
real a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
real b4=(k4m*pow((ksl/Kki),2.0))/(pow((1.0+nasl/Knai),3.0)+pow((1.0+ksl/Kki),2.0)-1.0);

x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;

E1=x1/(x1+x2+x3+x4);    
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);    
E4=x4/(x1+x2+x3+x4);

real JnakNa=3.0*(E1*a3-E2*b3);   
real JnakK=2.0*(E4*b1-E3*a1);    
real Pnak=32.4872;
Pnak *= GNaK_Multiplier;
real INaK=Pnak*(zna*JnakNa+zk*JnakK);
INaK *= INaK_Multiplier;

// Background currents: INab, ICab [from ORd]
// INab current
real PNab = 3.75e-10*2.5;
real INab = PNab*vffrt*(nasl*exp(vfrt)-nao)/(exp(vfrt)-1.0);
INab *= INab_Multiplier;
// ICab current
real PCab = 2.5e-8;
real ICab = PCab*4.0*vffrt*(casl*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
ICab *= ICab_Multiplier;
// IpCa current  [from PRd]
real ipcabar = 0.0005; 
real kmpca = 0.0005;
real IpCa= ipcabar/((kmpca/casl)+1);
IpCa *= IpCa_Multiplier; 
// Stimulation
real Istim = calc_I_stim;
//real dv = - (INa+INaL+Ito+Isus+ICaL+ICaNa+ICaK+ICaT+If+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IpCa+ICab+Istim);
real dv = - (INa+INaL+Ito+Isus+ICaL+ICaNa+ICaK+ICaT+If+IKr+IKs+IK1+INaCa+INaK+INab+IpCa+ICab+Istim);

// ------------------------------------------------------------------------
//                        SR_CA_FLUXES [from PRd]                      
// ------------------------------------------------------------------------
// Diffusion Fluxes
real JdiffNa= (nass-nasl)/sstau;
real JgapNa=(nasl-nai)/gaptau;
real JdiffK= (kss-ksl)/sstau;
real JgapK=(ksl-ki)/gaptau;
real Jdiff = (cass-casl)/sstau;  
real Jgap  = (casl-cai)/gaptau;  

// Ip3R_Ca_Release
real IP3 = 0.0001;  
k1 = 150000;
real k1a = 16.5;
real k0 = 96000;
real k0a = 9.6;
k2 = 1800;
real k2a = 0.21;
real tauip3r = 3.7;
real du=(cass*k2*(1-u)-k2a*u);
real POip3 = tauip3r*IP3*cass*(1-u)/((1+IP3*k0/k0a)*(1+cass*k1/k1a));
real Jip3 = 10.920*(cajsr-cass)*(POip3);                    

// Ca_uptake_via_SERCA: qup1 & qup2
real dqupcamkbar = 0.75;			        
real dkmplbbar = 0.00017;	
real kmup   = 0.00028;	
real nsrbar = 15.0;	 
real dkmplb		= dkmplbbar*CaMKa/(KmCaMK+CaMKa);
real dqupcamk	= dqupcamkbar*CaMKa/(KmCaMK+CaMKa); 
real Jup1		= (0.0002*(dqupcamk+1)/(1+((kmup-dkmplb)/casl))-0.00105*cansr/nsrbar);
real Jup2		= (0.0026*(dqupcamk+1)/(1+((kmup-dkmplb)/cai))-0.0042*cansr/nsrbar);

// RyR3_Ca_Release: qrel1 
real REL1  = -((ICaL)*Acap/(vss*2.0*F)-(Jrel1 + Jip3)*vjsr/vss + Jdiff);
real ireltau1 = 2*(1+1*(1/(1+(pow((0.28/CaMKa),8)))))/(1+(0.0123/cajsr));
real irelss = (REL1 > 0) ? 15*(1+1*(1/(1+(pow((0.28/CaMKa),8)))))*REL1/(1 + (pow((1.0/cajsr),8))) : 0;
real dJrel1= ((irelss-Jrel1)/ireltau1);

// RyR2_Ca_release: qrel2 
real REL2  = (-Jup2*vnsr/vmyo + Jgap*vsl/vmyo+ (Jrel2)*vcsr/vmyo);    
real ireltau = 6*(1+1*(1/(1+(pow((0.28/CaMKa),8)))))/(1+(0.0123/cacsr));
irelss = (REL2 > 0) ? 91*(1+1*(1/(1+(pow((0.28/CaMKa),8)))))*(REL2)/(1 + (pow((1/cacsr),8))) : 0;
real dJrel2= ((irelss-Jrel2)/ireltau);

// Tranlocation Flux
real Jtr1=(cansr-cajsr)/tautr1;
real Jtr2=(cansr-cacsr)/tautr2;

// Calcium Buffer Constants
//BSR
real BSRmax=0.019975;
real KmBSR=0.00087;
//BSL
real BSLmax=0.4777;
real KmBSL=0.0087;
//CSQN
real csqnmax= 2.88;    
real kmcsqn=0.8;     
real csqnmax1 = 1.2;       
//CMDN
real cmdnmax=0.1125;   
real kmcmdn=0.00238;
real cmdnmax1 = 1.25e-2;   
//TRPM
real trpnmax=3.15e-2;
real kmtrpn=0.0005;
real trpnmax1 = 3.5e-3;

// ------------------------------------------------------------------------
//                       IONIC CONCENTRATIONS                          
// ------------------------------------------------------------------------
// [Na]
real dnass=-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss*zna)-JdiffNa;
real dnasl=-(INa+INaL+3.0*INaCa_i+3.0*INaK+IfNa+INab)*Acap/(F*vsl*zna)+JdiffNa*vss/vsl-JgapNa;
real dnai=JgapNa*vsl/vmyo;
// [K]
real dkss=-(ICaK)*Acap/(F*vss*zk)-JdiffK;
real dksl=-(Ito+Isus+IKr+IKs+IfK+IK1+Istim-2.0*INaK)*Acap/(F*vsl*zk)+JdiffK*vss/vsl -JgapK;
real dki=JgapK*vsl/vmyo;
// Intracellular [Ca]
// i (Myo)
real trpnMYO=trpnmax*kmtrpn/pow((kmtrpn+cai),2.0);
real cmdnMYO=cmdnmax*kmcmdn/pow((kmcmdn+cai),2.0);
real Bcai   = 1.0 / (1.0+cmdnMYO+trpnMYO);
real dcaibar= (-(Jup2)*vnsr/vmyo+Jgap*vsl/vmyo+(Jrel2)*vcsr/vmyo);
real dcai = Bcai*dcaibar;
// SS
real Bcass  = 1.0/(1.0+BSRmax*KmBSR/pow((KmBSR+cass),2.0)+BSLmax*KmBSL/pow((KmBSL+cass),2.0));
real dcass  =  Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+(Jrel1+Jip3)*vjsr/vss-Jdiff);
// SL
real trpnSSL=(trpnmax1*kmtrpn)/(pow((kmtrpn+casl),2));
real cmdnSSL=(cmdnmax1*kmcmdn)/(pow((kmcmdn+casl),2));
real Bcasl=1/(1+trpnSSL+cmdnSSL);
real dcaslbar= (-(Jup1)*vnsr/vsl+Jdiff*vss/vsl-Jgap-(IpCa+ICab+ICaT-2*INaCa_i)*Acap/(vsl*2.0*F));
real dcasl=Bcasl*dcaslbar;
// Sarcoplasmic [Ca]
// NSR
real dcansr= Jup1+Jup2-Jtr1*vjsr/vnsr-Jtr2*vcsr/vnsr;
// JSR
real csqn1=(csqnmax1*kmcsqn)/(pow((kmcsqn+cajsr),2));
real Bcajsr = 1.0/(1.0+csqn1);
real dcajsrbar=(Jtr1-Jrel1-Jip3);
real dcajsr = Bcajsr*dcajsrbar;
// CSR	
real csqn=(csqnmax*kmcsqn)/(pow((cacsr+kmcsqn),2));
real Bcacsr=1/(1+csqn);
real dcacsrbar=(Jtr2-Jrel2);
real dcacsr=Bcacsr*dcacsrbar;

rDY_[0] = dv;      
rDY_[1] = dCaMKt;  
rDY_[2] = dcass;
rDY_[3] = dnai;    
rDY_[4] = dnasl;   
rDY_[5] = dnass;
rDY_[6] = dki; 
rDY_[8] = dksl;    
rDY_[7] = dkss;    
rDY_[9] = dcai;
rDY_[10] = dcasl;   
rDY_[11] = dcansr;
rDY_[12] = dcajsr;
rDY_[13] = dcacsr;  
rDY_[14] = dJrel1;
rDY_[15] = dJrel2;
rDY_[16] = dm;      
rDY_[17] = dhf;
rDY_[18] = dhs;     
rDY_[19] = dj;      
rDY_[20] = dhsp;
rDY_[21] = djp;     
rDY_[22] = dmL;
rDY_[23] = dhL;
rDY_[24] = dhLp;    
rDY_[25] = da;      
rDY_[26] = di1;   
rDY_[27] = di2;
rDY_[28] = dd;      
rDY_[29] = dff;     
rDY_[30] = dfs;
rDY_[31] = dfcaf;   
rDY_[32] = dfcas;
rDY_[33] = djca;    
rDY_[34] = dffp;    
rDY_[35] = dfcafp;  
rDY_[36] = dnca;    
rDY_[37] = db;
rDY_[38] = dg;      
rDY_[39] = dxrf;    
rDY_[40] = dxrs;    
rDY_[41] = dxs1;    
rDY_[42] = dxs2;
rDY_[43] = dy;      
rDY_[44] = dxk1;    
rDY_[45] = du;      