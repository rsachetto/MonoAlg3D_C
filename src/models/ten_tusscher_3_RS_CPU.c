#include "model_common.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#define ENDO

void RHS_cpu(const Real *sv, Real *rDY_, Real stim_current, Real time, Real dt);


void init_cell_model_data(struct cell_model_data* cell_model, bool get_initial_v, bool get_neq) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = -86.2f;
    if(get_neq)
        cell_model->number_of_ode_equations = 12;

}

void set_model_initial_conditions_cpu(Real *sv) {

    sv[0] = -85.23f;   // V;       millivolt
    sv[1] = 0.; //M
    sv[2] = 0.75; //H
    sv[3] = 0.75; //J
    sv[4] = 0.; //Xr1
    sv[5] = 0.; //Xs
    sv[6] = 1.; //S
    sv[7] = 1.; //F
    sv[8] = 1.; //F2
    sv[9] = 0.0; //D_INF
    sv[10] = 0.0; //R_INF
    sv[11] = 0.0; //Xr2_INF}
}

void solve_model_ode_cpu(Real dt, Real *sv, Real stim_current, Real time, int neq, void *extra_data)  {

    assert(sv);
    extra_data = NULL;

    Real rY[neq], rDY[neq];

    for(int i = 0; i < neq; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, time, dt);

    //THIS MODEL USES THE Rush Larsen Method TO SOLVE THE EDOS
    sv[0] = dt*rDY[0] + rY[0];
    sv[1]  = rDY[1];
    sv[2]  = rDY[2];
    sv[3]  = rDY[3];
    sv[4]  = rDY[4];
    sv[5]  = rDY[5];
    sv[6]  = rDY[6];
    sv[7]  = rDY[7];
    sv[8]  = rDY[8];
    sv[9]  = rDY[9];
    sv[10]  = rDY[10];
    sv[11]  = rDY[11];
}


void RHS_cpu(const Real *sv, Real *rDY_, Real stim_current, Real time, Real dt) {

    // State variables
    const Real svolt = sv[0];      // Membrane variable

    const Real sm   = sv[1];
    const Real sh   = sv[2];
    const Real sj   = sv[3];
    const Real sxr1 = sv[4];
    const Real sxs  = sv[5];
    const Real ss   = sv[6];
    const Real  sf   = sv[7];
    const Real sf2  = sv[8];

    const Real D_INF  = sv[9];
    const Real Xr2_INF  = sv[10];
    const Real R_INF  = sv[11];


    const Real natp = 0.24;          // K dependence of ATP-sensitive K current
    const Real nicholsarea = 0.00005; // Nichol's areas (cm^2)
    //const Real atpi = 4.0;             // Intracellular ATP concentraion (mM)
    const Real hatp = 2;             // Hill coefficient
    //const Real katp = 0.306;         // Half-maximal saturation point of ATP-sensitive K current (mM)


    const Real  atpi = Atpi;
    const Real katp = KAtp;

    const Real patp =  1/(1 + pow((atpi/katp),hatp));
    const Real gkatp    =  0.000195/nicholsarea;
    const Real gkbaratp =  gkatp*patp*pow((Ko/4),natp);


    const Real Cao=2.0;
    const Real Nao=140.0;
    const Real Cai=0.00007;
    const Real Nai=7.67;
    const Real Ki=138.3;

//Constants
    const Real R=8314.472;
    const Real F=96485.3415;
    const Real T=310.0;
    const Real RTONF=(R*T)/F;

//Parameters for currents
//Parameters for IKr
    const Real Gkr=0.101;
//Parameters for Iks
    const Real pKNa=0.03;
#ifdef EPI
    const Real Gks=0.257;
#endif
#ifdef ENDO
    const Real Gks=0.392;
#endif
#ifdef MCELL
    const Real Gks=0.098;
#endif
//Parameters for Ik1
    const Real GK1=5.405;
//Parameters for Ito
#ifdef EPI
    const Real Gto=0.294;
#endif
#ifdef ENDO
    const Real Gto=0.073;
#endif
#ifdef MCELL
    const Real Gto=0.294;
#endif
//Parameters for INa
    const Real GNa=14.838;
//Parameters for IbNa
    const Real GbNa=0.00029;
//Parameters for INaK
    const Real KmK=1.0;
    const Real KmNa=40.0;
    const Real knak=2.724;
//Parameters for ICaL
    const Real GCaL=0.2786;
//Parameters for IbCa
    const Real GbCa=0.000592;
//Parameters for INaCa
    const Real knaca=1000;
    const Real KmNai=87.5;
    const Real KmCa=1.38;
    const Real ksat=0.1;
    const Real n=0.35;
//Parameters for IpCa
    const Real GpCa=0.1238;
    const Real KpCa=0.0005;
//Parameters for IpK;
    const Real GpK=0.0293;


    Real Ek=RTONF*(log((Ko/Ki)));
    Real Ena=RTONF*(log((Nao/Nai)));
    Real Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    Real Eca=0.5*RTONF*(log((Cao/Cai)));
    Real IKr;
    Real IKs;
    Real IK1;
    Real Ito;
    Real INa;
    Real IbNa;
    Real ICaL;
    Real IbCa;
    Real INaCa;
    Real IpCa;
    Real IpK;
    Real INaK;
    Real IKatp;

    static Real Ak1;
    static Real Bk1;
    static Real rec_iK1;
    static Real rec_ipK;
    static Real rec_iNaK;
    static Real AM;
    static Real BM;
    static Real AH_1;
    static Real BH_1;
    static Real AH_2;
    static Real BH_2;
    static Real AJ_1;
    static Real BJ_1;
    static Real AJ_2;
    static Real BJ_2;
    static Real M_INF;
    static Real H_INF;
    static Real J_INF;
    static Real TAU_M;
    static Real TAU_H;
    static Real TAU_J;
    static Real axr1;
    static Real bxr1;
    static Real Xr1_INF;
    static Real Xr2_INF_new;
    static Real TAU_Xr1;
    static Real Axs;
    static Real Bxs;
    static Real Xs_INF;
    static Real TAU_Xs;
    static Real R_INF_new;
    static Real S_INF;
    static Real TAU_S;
    static Real Af;
    static Real Bf;
    static Real Cf;
    static Real Af2;
    static Real Bf2;
    static Real Cf2;
    static Real D_INF_new;
    static Real TAU_F;
    static Real F_INF;
    static Real TAU_F2;
    static Real F2_INF;
    Real sItot;


    //Needed to compute currents
    Ak1=0.1/(1.+exp(0.06*(svolt-Ek-200)));
    Bk1=(3.*exp(0.0002*(svolt-Ek+100))+
         exp(0.1*(svolt-Ek-10)))/(1.+exp(-0.5*(svolt-Ek)));
    rec_iK1=Ak1/(Ak1+Bk1);
    rec_iNaK=(1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
    rec_ipK=1./(1.+exp((25-svolt)/5.98));


    //Compute currents
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
    ICaL=GCaL*D_INF*sf*sf2*(svolt-60);
    Ito=Gto*R_INF*ss*(svolt-Ek);
    IKr=Gkr*sqrt(Ko/5.4)*sxr1*Xr2_INF*(svolt-Ek);
    IKs=Gks*sxs*sxs*(svolt-Eks);
    IK1=GK1*rec_iK1*(svolt-Ek);
    INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
          (1./(1+ksat*exp((n-1)*svolt*F/(R*T))))*
          (exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
           exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa=GpCa*Cai/(KpCa+Cai);
    IpK=GpK*rec_ipK*(svolt-Ek);
    IbNa=GbNa*(svolt-Ena);
    IbCa=GbCa*(svolt-Eca);

    IKatp = gkbaratp*(svolt-Ek);

    Real calc_i_stim = ((time>=stim_start)&&(time<=stim_start+stim_dur)) ? Istim: 0.0;


    //Determine total current
    (sItot) = IKr    +
              IKs   +
              IK1   +
              Ito   +
              INa   +
              IbNa  +
              ICaL  +
              IbCa  +
              INaK  +
              INaCa +
              IpCa  +
              IpK   +
              IKatp +
              calc_i_stim;

    //compute steady state values and time constants
    AM=1./(1.+exp((-60.-svolt)/5.));
    BM=0.1/(1.+exp((svolt+35.)/5.))+0.10/(1.+exp((svolt-50.)/200.));
    TAU_M=AM*BM;
    M_INF=1./((1.+exp((-56.86-svolt)/9.03))*(1.+exp((-56.86-svolt)/9.03)));
    if (svolt>=-40.)
    {
        AH_1=0.;
        BH_1=(0.77/(0.13*(1.+exp(-(svolt+10.66)/11.1))));
        TAU_H= 1.0/(AH_1+BH_1);
    }
    else
    {
        AH_2=(0.057*exp(-(svolt+80.)/6.8));
        BH_2=(2.7*exp(0.079*svolt)+(3.1e5)*exp(0.3485*svolt));
        TAU_H=1.0/(AH_2+BH_2);
    }
    H_INF=1./((1.+exp((svolt+71.55)/7.43))*(1.+exp((svolt+71.55)/7.43)));
    if(svolt>=-40.)
    {
        AJ_1=0.;
        BJ_1=(0.6*exp((0.057)*svolt)/(1.+exp(-0.1*(svolt+32.))));
        TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
        AJ_2=(((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)*
                                             exp(-0.04391*svolt))*(svolt+37.78)/
              (1.+exp(0.311*(svolt+79.23))));
        BJ_2=(0.02424*exp(-0.01052*svolt)/(1.+exp(-0.1378*(svolt+40.14))));
        TAU_J= 1.0/(AJ_2+BJ_2);
    }
    J_INF=H_INF;

    Xr1_INF=1./(1.+exp((-26.-svolt)/7.));
    axr1=450./(1.+exp((-45.-svolt)/10.));
    bxr1=6./(1.+exp((svolt-(-30.))/11.5));
    TAU_Xr1=axr1*bxr1;
    Xr2_INF_new=1./(1.+exp((svolt-(-88.))/24.));


    Xs_INF=1./(1.+exp((-5.-svolt)/14.));
    Axs=(1400./(sqrt(1.+exp((5.-svolt)/6))));
    Bxs=(1./(1.+exp((svolt-35.)/15.)));
    TAU_Xs=Axs*Bxs+80;

#ifdef EPI
    R_INF_new=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+20)/5.));
    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
#endif
#ifdef ENDO
    R_INF_new=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+28)/5.));
    TAU_S=1000.*exp(-(svolt+67)*(svolt+67)/1000.)+8.;
#endif
#ifdef MCELL
    R_INF_new=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+20)/5.));
    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
#endif


    D_INF_new=1./(1.+exp((-8-svolt)/7.5));
    F_INF=1./(1.+exp((svolt+20)/7));
    Af=1102.5*exp(-(svolt+27)*(svolt+27)/225);
    Bf=200./(1+exp((13-svolt)/10.));
    Cf=(180./(1+exp((svolt+30)/10)))+20;
    TAU_F=Af+Bf+Cf;
    F2_INF=0.67/(1.+exp((svolt+35)/7))+0.33;
    Af2=600*exp(-(svolt+27)*(svolt+27)/170);
    Bf2=7.75/(1.+exp((25-svolt)/10));
    Cf2=16/(1.+exp((svolt+30)/10));
    TAU_F2=Af2+Bf2+Cf2;

    //update voltage
    rDY_[0] = -sItot;

    //Update gates
    const Real dt = global_dt;
    rDY_[1] = M_INF-(M_INF-sm)*exp(-dt/TAU_M);
    rDY_[2] = H_INF-(H_INF-sh)*exp(-dt/TAU_H);
    rDY_[3] = J_INF-(J_INF-sj)*exp(-dt/TAU_J);
    rDY_[4] = Xr1_INF-(Xr1_INF-sxr1)*exp(-dt/TAU_Xr1);
    rDY_[5] = Xs_INF-(Xs_INF-sxs)*exp(-dt/TAU_Xs);
    rDY_[6]= S_INF-(S_INF-ss)*exp(-dt/TAU_S);
    rDY_[7] =F_INF-(F_INF-sf)*exp(-dt/TAU_F);
    rDY_[8] =F2_INF-(F2_INF-sf2)*exp(-dt/TAU_F2);

    rDY_[9] = D_INF_new;
    rDY_[10] = R_INF_new;
    rDY_[11] = Xr2_INF_new;

}