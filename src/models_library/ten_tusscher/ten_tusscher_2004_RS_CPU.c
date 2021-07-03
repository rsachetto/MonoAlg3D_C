#include <assert.h>
#include <stdlib.h>
#include "ten_tusscher_2004.h"


GET_CELL_MODEL_DATA(init_cell_model_data) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;

}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {

        real *sv = &solver->sv[i * NEQ];

        sv[0] = INITIAL_V; // V;       millivolt
        sv[1] = 0.f;       // M
        sv[2] = 0.75;      // H
        sv[3] = 0.75f;     // J
        sv[4] = 0.f;       // Xr1
        sv[5] = 1.f;       // Xr2
        sv[6] = 0.f;       // Xs
        sv[7] = 1.f;       // S
        sv[8] = 0.f;       // R
        sv[9] = 0.f;       // D
        sv[10] = 1.f;      // F
        sv[11] = 1.f;      // FCa
        sv[12] = 1.f;      // G
        sv[13] = 0.0002;   // Cai
        sv[14] = 0.2f;     // CaSR
        sv[15] = 11.6f;    // Nai
        sv[16] = 138.3f;   // Ki
    }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;
    int i;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;

    OMP(parallel for private(sv_id))
    for (i = 0; i < num_cells_to_solve; i++) {

        if(ode_solver->cells_to_solve)
            sv_id = ode_solver->cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < ode_solver->num_steps; ++j) {
            solve_model_ode_cpu(ode_solver->min_dt, ode_solver->sv + (sv_id * NEQ), stim_currents[i]);

        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt);

    for(int i = 0; i < NEQ; i++)
        sv[i] = rDY[i];
}


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {


    // State variables
    real svolt = sv[0];
    real sm    = sv[1];
    real sh    = sv[2];
    real sj    = sv[3];
    real sxr1  = sv[4];
    real sxr2  = sv[5];
    real sxs   = sv[6];
    real ss    = sv[7];
    real sr    = sv[8];
    real sd    = sv[9];
    real sf    = sv[10];
    real sfca  = sv[11];
    real sg    = sv[12];
    real Cai   = sv[13];
    real CaSR  = sv[14];
    real Nai   = sv[15];
    real Ki    = sv[16];

    //External concentrations
    real Ko=5.4;
    real Cao=2.0;
    real Nao=140.0;

    //Intracellular volumes
    real Vc=0.016404;
    real Vsr=0.001094;

    //Calcium dynamics
    real Bufc=0.15f;
    real Kbufc=0.001f;
    real Bufsr=10.f;
    real Kbufsr=0.3f;
    real taufca=2.f;
    real taug=2.f;
    real Vmaxup=0.000425f;
    real Kup=0.00025f;

//Constants
    const real R = 8314.472f;
    const real F = 96485.3415f;
    const real T =310.0f;
    real RTONF   =(R*T)/F;

//Cellular capacitance
    real CAPACITANCE=0.185;

//Parameters for currents
//Parameters for IKr
    real Gkr=0.096;
//Parameters for Iks
    real pKNa=0.03;
#ifdef EPI
    real Gks=0.245;
#endif
#ifdef ENDO
    real Gks=0.245;
#endif
#ifdef MCELL
    real Gks=0.062;
#endif
//Parameters for Ik1
    real GK1=5.405;
//Parameters for Ito
#ifdef EPI
    real Gto=0.294;
#endif
#ifdef ENDO
    real Gto=0.073;
#endif
#ifdef MCELL
    real Gto=0.294;
#endif
//Parameters for INa
    real GNa=14.838;
//Parameters for IbNa
    real GbNa=0.00029;
//Parameters for INaK
    real KmK=1.0;
    real KmNa=40.0;
    real knak=1.362;
//Parameters for ICaL
    real GCaL=0.000175;
//Parameters for IbCa
    real GbCa=0.000592;
//Parameters for INaCa
    real knaca=1000;
    real KmNai=87.5;
    real KmCa=1.38;
    real ksat=0.1;
    real n=0.35;
//Parameters for IpCa
    real GpCa=0.825;
    real KpCa=0.0005;
//Parameters for IpK;
    real GpK=0.0146;


    real IKr;
    real IKs;
    real IK1;
    real Ito;
    real INa;
    real IbNa;
    real ICaL;
    real IbCa;
    real INaCa;
    real IpCa;
    real IpK;
    real INaK;
    real Irel;
    real Ileak;


    real dNai;
    real dKi;
    real dCai;
    real dCaSR;

    real A;
//    real BufferFactorc;
//    real BufferFactorsr;
    real SERCA;
    real Caisquare;
    real CaSRsquare;
    real CaCurrent;
    real CaSRCurrent;


    real fcaold;
    real gold;
    real Ek;
    real Ena;
    real Eks;
    real Eca;
    real CaCSQN;
    real bjsr;
    real cjsr;
    real CaBuf;
    real bc;
    real cc;
    real Ak1;
    real Bk1;
    real rec_iK1;
    real rec_ipK;
    real rec_iNaK;
    real AM;
    real BM;
    real AH_1;
    real BH_1;
    real AH_2;
    real BH_2;
    real AJ_1;
    real BJ_1;
    real AJ_2;
    real BJ_2;
    real M_INF;
    real H_INF;
    real J_INF;
    real TAU_M;
    real TAU_H;
    real TAU_J;
    real axr1;
    real bxr1;
    real axr2;
    real bxr2;
    real Xr1_INF;
    real Xr2_INF;
    real TAU_Xr1;
    real TAU_Xr2;
    real Axs;
    real Bxs;
    real Xs_INF;
    real TAU_Xs;
    real R_INF;
    real TAU_R;
    real S_INF;
    real TAU_S;
    real Ad;
    real Bd;
    real Cd;
    real TAU_D;
    real D_INF;
    real TAU_F;
    real F_INF;
    real FCa_INF;
    real G_INF;

    real inverseVcF2=1/(2*Vc*F);
    real inverseVcF=1./(Vc*F);
    real Kupsquare=Kup*Kup;
//    real BufcKbufc=Bufc*Kbufc;
//    real Kbufcsquare=Kbufc*Kbufc;
//    real Kbufc2=2*Kbufc;
//    real BufsrKbufsr=Bufsr*Kbufsr;
//    const real Kbufsrsquare=Kbufsr*Kbufsr;
//    const real Kbufsr2=2*Kbufsr;
    const real exptaufca=exp(-dt/taufca);
    const real exptaug=exp(-dt/taug);

    real sItot;

    //Needed to compute currents
    Ek=RTONF*(log((Ko/Ki)));
    Ena=RTONF*(log((Nao/Nai)));
    Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    Eca=0.5*RTONF*(log((Cao/Cai)));
    Ak1=0.1/(1.+exp(0.06*(svolt-Ek-200)));
    Bk1=(3.*exp(0.0002*(svolt-Ek+100))+
         exp(0.1*(svolt-Ek-10)))/(1.+exp(-0.5*(svolt-Ek)));
    rec_iK1=Ak1/(Ak1+Bk1);
    rec_iNaK=(1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
    rec_ipK=1./(1.+exp((25-svolt)/5.98));


    //Compute currents
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
    ICaL=GCaL*sd*sf*sfca*4*svolt*(F*F/(R*T))*
         (exp(2*svolt*F/(R*T))*Cai-0.341*Cao)/(exp(2*svolt*F/(R*T))-1.);
    Ito=Gto*sr*ss*(svolt-Ek);
    IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);
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
              stim_current;


    //update concentrations
    Caisquare=Cai*Cai;
    CaSRsquare=CaSR*CaSR;
    CaCurrent=-(ICaL+IbCa+IpCa-2.0f*INaCa)*inverseVcF2*CAPACITANCE;
    A=0.016464f*CaSRsquare/(0.0625f+CaSRsquare)+0.008232f;
    Irel=A*sd*sg;
    Ileak=0.00008f*(CaSR-Cai);
    SERCA=Vmaxup/(1.f+(Kupsquare/Caisquare));
    CaSRCurrent=SERCA-Irel-Ileak;
    CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    dCaSR=dt*(Vc/Vsr)*CaSRCurrent;
    bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
    cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
    CaSR=(sqrt(bjsr*bjsr+4.*cjsr)-bjsr)/2.;
    CaBuf=Bufc*Cai/(Cai+Kbufc);
    dCai=dt*(CaCurrent-CaSRCurrent);
    bc=Bufc-CaBuf-dCai-Cai+Kbufc;
    cc=Kbufc*(CaBuf+dCai+Cai);
    Cai=(sqrt(bc*bc+4*cc)-bc)/2;



    dNai=-(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
    Nai+=dt*dNai;

    dKi=-(stim_current+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
    Ki+=dt*dKi;

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
    Xr2_INF=1./(1.+exp((svolt-(-88.))/24.));
    axr2=3./(1.+exp((-60.-svolt)/20.));
    bxr2=1.12/(1.+exp((svolt-60.)/20.));
    TAU_Xr2=axr2*bxr2;

    Xs_INF=1./(1.+exp((-5.-svolt)/14.));
    Axs=1100./(sqrt(1.+exp((-10.-svolt)/6)));
    Bxs=1./(1.+exp((svolt-60.)/20.));
    TAU_Xs=Axs*Bxs;

#ifdef EPI
    R_INF=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+20)/5.));
    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
#endif
#ifdef ENDO
    R_INF=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+28)/5.));
    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=1000.*exp(-(svolt+67)*(svolt+67)/1000.)+8.;
#endif
#ifdef MCELL
    R_INF=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+20)/5.));
    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
#endif


    D_INF=1./(1.+exp((-5-svolt)/7.5));
    Ad=1.4/(1.+exp((-35-svolt)/13))+0.25;
    Bd=1.4/(1.+exp((svolt+5)/5));
    Cd=1./(1.+exp((50-svolt)/20));
    TAU_D=Ad*Bd+Cd;
    F_INF=1./(1.+exp((svolt+20)/7));
    TAU_F=1125*exp(-(svolt+27)*(svolt+27)/300)+80+165/(1.+exp((25-svolt)/10));


    FCa_INF=(1./(1.+pow((Cai/0.000325),8))+
             0.1/(1.+exp((Cai-0.0005)/0.0001))+
             0.20/(1.+exp((Cai-0.00075)/0.0008))+
             0.23 )/1.46;
    if(Cai<0.00035)
        G_INF=1./(1.+pow((Cai/0.00035),6));
    else
        G_INF=1./(1.+pow((Cai/0.00035),16));

    //Update gates
    rDY_[1]  = M_INF-(M_INF-sm)*exp(-dt/TAU_M);
    rDY_[2]  = H_INF-(H_INF-sh)*exp(-dt/TAU_H);
    rDY_[3]  = J_INF-(J_INF-sj)*exp(-dt/TAU_J);
    rDY_[4]  = Xr1_INF-(Xr1_INF-sxr1)*exp(-dt/TAU_Xr1);
    rDY_[5]  = Xr2_INF-(Xr2_INF-sxr2)*exp(-dt/TAU_Xr2);
    rDY_[6]  = Xs_INF-(Xs_INF-sxs)*exp(-dt/TAU_Xs);
    rDY_[7]  = S_INF-(S_INF-ss)*exp(-dt/TAU_S);
    rDY_[8]  = R_INF-(R_INF-sr)*exp(-dt/TAU_R);
    rDY_[9]  = D_INF-(D_INF-sd)*exp(-dt/TAU_D);
    rDY_[10] = F_INF-(F_INF-sf)*exp(-dt/TAU_F);
    fcaold= sfca;
    sfca = FCa_INF-(FCa_INF-sfca)*exptaufca;
    if(sfca>fcaold && (svolt)>-37.0)
        sfca = fcaold;
    gold = sg;
    sg = G_INF-(G_INF-sg)*exptaug;

    if(sg>gold && (svolt)>-37.0)
        sg=gold;

    //update voltage
    rDY_[0] = svolt + dt*(-sItot);
    rDY_[11] = sfca;
    rDY_[12] = sg;
    rDY_[13] = Cai;
    rDY_[14] = CaSR;
    rDY_[15] = Nai;
    rDY_[16] = Ki;

}
