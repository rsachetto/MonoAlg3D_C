#include "model_common.h"
#include <assert.h>
#include <stdlib.h>
#include "ten_tusscher_3_RS.h"


GET_CELL_MODEL_DATA(init_cell_model_data) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;

}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    sv[0] = -86.2f;   // V;       millivolt
    sv[1] = 0.0f; //M
    sv[2] = 0.75; //H
    sv[3] = 0.75; //J
    sv[4] = 0.0f; //Xr1
    sv[5] = 0.0f; //Xs
    sv[6] = 1.0f; //S
    sv[7] = 1.0f; //F
    sv[8] = 1.0f; //F2
    sv[9] = 0.0; //D_INF
    sv[10] = 0.0; //R_INF
    sv[11] = 0.0; //Xr2_INF}
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    uint32_t sv_id;
    real atpi;
    real *fibrosis;

    if(extra_data) {
        atpi = *((real*)extra_data);
        fibrosis = ((real*)extra_data)+1;
    }
    else {
        atpi = 6.8;
        fibrosis=calloc(num_cells_to_solve, sizeof(real));
    }


	int i;

    #pragma omp parallel for private(sv_id)
    for (i = 0; i < num_cells_to_solve; i++) {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], fibrosis[i], atpi);

        }
    }

    if(extra_data == NULL) free(fibrosis);
}


void solve_model_ode_cpu(real dt, real *sv, real stim_current, real fibrosis, real atpi )  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, fibrosis, atpi);

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


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real fibrosis, real atpi) {

    // State variables
    const real svolt = sv[0];      // Membrane variable

    const real sm   = sv[1];
    const real sh   = sv[2];
    const real sj   = sv[3];
    const real sxr1 = sv[4];
    const real sxs  = sv[5];
    const real ss   = sv[6];
    const real  sf   = sv[7];
    const real sf2  = sv[8];

    const real D_INF  = sv[9];
    const real Xr2_INF  = sv[10];
    const real R_INF  = sv[11];


    const real natp = 0.24;          // K dependence of ATP-sensitive K current
    const real nicholsarea = 0.00005; // Nichol's areas (cm^2)
    const real hatp = 2;             // Hill coefficient

    real Ko   = 5.4;

    real atpi_change = 6.8f - atpi;

    atpi = atpi + atpi_change*fibrosis;

    //real katp = 0.306;
    const real katp = -0.0942857142857f*atpi + 0.683142857143f; //Ref: A Comparison of Two Models of Human Ventricular Tissue: Simulated Ischaemia and Re-entry


    const real patp =  1.0f/(1.0f + powf((atpi/katp),hatp));
    const real gkatp    =  0.000195f/nicholsarea;
    const real gkbaratp =  gkatp*patp*powf((Ko/4),natp);

    const real katp2= 1.4;
    const real hatp2 = 2.6;
    const real pcal = 1.0f/(1.0f + powf((katp2/atpi),hatp2));

    const real Cao=2.0;
    const real Nao=140.0;
    const real Cai=0.00007;
    const real Nai=7.67;
    const real Ki=138.3;

//Constants
    const real R=8314.472;
    const real F=96485.3415f;
    const real T=310.0;
    const real RTONF=(R*T)/F;

//Parameters for currents
//Parameters for IKr
    const real Gkr=0.101;
//Parameters for Iks
    const real pKNa=0.03;
#ifdef EPI
    const real Gks=0.257;
#endif
#ifdef ENDO
    const real Gks=0.392;
#endif
#ifdef MCELL
    const real Gks=0.098;
#endif
//Parameters for Ik1
    const real GK1=5.405;
//Parameters for Ito
#ifdef EPI
    const real Gto=0.294;
#endif
#ifdef ENDO
    const real Gto=0.073;
#endif
#ifdef MCELL
    const real Gto=0.294;
#endif
//Parameters for INa
    const real GNa=14.838;
//Parameters for IbNa
    const real GbNa=0.00029;
//Parameters for INaK
    const real KmK=1.0;
    const real KmNa=40.0;
    const real knak=2.724;
//Parameters for ICaL
    const real GCaL=0.2786f*pcal;
//Parameters for IbCa
    const real GbCa=0.000592;
//Parameters for INaCa
    const real knaca=1000;
    const real KmNai=87.5;
    const real KmCa=1.38;
    const real ksat=0.1;
    const real n=0.35;
//Parameters for IpCa
    const real GpCa=0.1238;
    const real KpCa=0.0005;
//Parameters for IpK;
    const real GpK=0.0293;


    const real Ek=RTONF*(logf((Ko/Ki)));
    const real Ena=RTONF*(logf((Nao/Nai)));
    const real Eks=RTONF*(logf((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    const real Eca=0.5f*RTONF*(logf((Cao/Cai)));
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
    real IKatp;

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
    real Xr1_INF;
    real Xr2_INF_new;
    real TAU_Xr1;
    real Axs;
    real Bxs;
    real Xs_INF;
    real TAU_Xs;
    real R_INF_new;
    real S_INF;
    real TAU_S;
    real Af;
    real Bf;
    real Cf;
    real Af2;
    real Bf2;
    real Cf2;
    real D_INF_new;
    real TAU_F;
    real F_INF;
    real TAU_F2;
    real F2_INF;
    real sItot;


    //Needed to compute currents
    Ak1=0.1f/(1.0f+expf(0.06f*(svolt-Ek-200.0f)));
    Bk1=(3.0f*expf(0.0002f*(svolt-Ek+100.0f))+
         expf(0.1f*(svolt-Ek-10.0f)))/(1.0f+expf(-0.5f*(svolt-Ek)));
    rec_iK1=Ak1/(Ak1+Bk1);
    rec_iNaK=(1.0f/(1.0f+0.1245f*expf(-0.1f*svolt*F/(R*T))+0.0353f*expf(-svolt*F/(R*T))));
    rec_ipK=1.0f/(1.0f+expf((25.0f-svolt)/5.98f));


    //Compute currents
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
    ICaL=GCaL*D_INF*sf*sf2*(svolt-60);
    Ito=Gto*R_INF*ss*(svolt-Ek);
    IKr=Gkr*sqrtf(Ko/5.4f)*sxr1*Xr2_INF*(svolt-Ek);
    IKs=Gks*sxs*sxs*(svolt-Eks);
    IK1=GK1*rec_iK1*(svolt-Ek);
    INaCa=knaca*(1.0f/(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1.0f/(KmCa+Cao))*
          (1.0f/(1.0f+ksat*expf((n-1.0f)*svolt*F/(R*T))))*
          (expf(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
           expf((n-1.0f)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5f);
    INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa=GpCa*Cai/(KpCa+Cai);
    IpK=GpK*rec_ipK*(svolt-Ek);
    IbNa=GbNa*(svolt-Ena);
    IbCa=GbCa*(svolt-Eca);

    IKatp = gkbaratp*(svolt-Ek);

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
              stim_current;

    //compute steady state values and time constants
    AM=1.0f/(1.0f+expf((-60.0f-svolt)/5.0f));
    BM=0.1f/(1.0f+expf((svolt+35.0f)/5.0f))+0.10f/(1.0f+expf((svolt-50.0f)/200.0f));
    TAU_M=AM*BM;
    M_INF=1.0f/((1.0f+expf((-56.86f-svolt)/9.03f))*(1.0f+expf((-56.86f-svolt)/9.03f)));
    if (svolt>=-40.)
    {
        AH_1=0.0f;
        BH_1=(0.77f/(0.13f*(1.0f+expf(-(svolt+10.66f)/11.1f))));
        TAU_H= 1.0f/(AH_1+BH_1);
    }
    else
    {
        AH_2=(0.057f*expf(-(svolt+80.0f)/6.8f));
        BH_2=(2.7f*expf(0.079f*svolt)+(3.1e5f)*expf(0.3485f*svolt));
        TAU_H=1.0f/(AH_2+BH_2);
    }
    H_INF=1.0f/((1.0f+expf((svolt+71.55f)/7.43f))*(1.0f+expf((svolt+71.55f)/7.43f)));
    if(svolt>=-40.0f)
    {
        AJ_1=0.0f;
        BJ_1=(0.6f*expf((0.057f)*svolt)/(1.0f+expf(-0.1f*(svolt+32.0f))));
        TAU_J= 1.0f/(AJ_1+BJ_1);
    }
    else
    {
        AJ_2=(((-2.5428e4f)*expf(0.2444f*svolt)-(6.948e-6f)*expf(-0.04391f*svolt))*(svolt+37.78f)/
              (1.0f+expf(0.311f*(svolt+79.23f))));
        BJ_2=(0.02424f*expf(-0.01052f*svolt)/(1.0f+expf(-0.1378f*(svolt+40.14f))));
        TAU_J= 1.0f/(AJ_2+BJ_2);
    }
    J_INF=H_INF;

    Xr1_INF=1.0f/(1.0f+expf((-26.0f-svolt)/7.0f));
    axr1=450.0f/(1.0f+expf((-45.0f-svolt)/10.0f));
    bxr1=6.0f/(1.0f+expf((svolt-(-30.0f))/11.5f));
    TAU_Xr1=axr1*bxr1;
    Xr2_INF_new=1.0f/(1.0f+expf((svolt-(-88.0f))/24.0f));


    Xs_INF=1.0f/(1.0f+expf((-5.0f-svolt)/14.0f));
    Axs=(1400.0f/(sqrtf(1.0f+expf((5.0f-svolt)/6.0f))));
    Bxs=(1.0f/(1.0f+expf((svolt-35.0f)/15.0f)));
    TAU_Xs=Axs*Bxs+80;

#ifdef EPI
    R_INF_new=1./(1.+expf((20-svolt)/6.));
    S_INF=1./(1.+expf((svolt+20)/5.));
    TAU_S=85.*expf(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+expf((svolt-20.)/5.))+3.;
#endif
#ifdef ENDO
    R_INF_new=1.0f/(1.0f+expf((20.0f-svolt)/6.0f));
    S_INF=1.0f/(1.0f+expf((svolt+28.0f)/5.0f));
    TAU_S=1000.0f*expf(-(svolt+67.0f)*(svolt+67.0f)/1000.0f)+8.0f;
#endif
#ifdef MCELL
    R_INF_new=1./(1.+expf((20-svolt)/6.));
    S_INF=1./(1.+expf((svolt+20)/5.));
    TAU_S=85.*expf(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+expf((svolt-20.)/5.))+3.;
#endif


    D_INF_new=1.0f/(1.0f+expf((-8.0f-svolt)/7.5f));
    F_INF=1.0f/(1.0f+expf((svolt+20)/7));
    Af=1102.5f*expf(-(svolt+27)*(svolt+27.0f)/225.0f);
    Bf=200.0f/(1.0f+expf((13.0f-svolt)/10.f));
    Cf=(180.0f/(1.0f+expf((svolt+30.0f)/10.0f)))+20.0f;
    TAU_F=Af+Bf+Cf;
    F2_INF=0.67f/(1.0f+expf((svolt+35.0f)/7.0f))+0.33f;
    Af2=600.0f*expf(-(svolt+27.0f)*(svolt+27.0f)/170.0f);
    Bf2=7.75f/(1.0f+expf((25.0f-svolt)/10.0f));
    Cf2=16.0f/(1.0f+expf((svolt+30.0f)/10.0f));
    TAU_F2=Af2+Bf2+Cf2;

    //update voltage
    rDY_[0] = -sItot;

    //Update gates
    rDY_[1] = M_INF-(M_INF-sm)*expf(-dt/TAU_M);
    rDY_[2] = H_INF-(H_INF-sh)*expf(-dt/TAU_H);
    rDY_[3] = J_INF-(J_INF-sj)*expf(-dt/TAU_J);
    rDY_[4] = Xr1_INF-(Xr1_INF-sxr1)*expf(-dt/TAU_Xr1);
    rDY_[5] = Xs_INF-(Xs_INF-sxs)*expf(-dt/TAU_Xs);
    rDY_[6]= S_INF-(S_INF-ss)*expf(-dt/TAU_S);
    rDY_[7] =F_INF-(F_INF-sf)*expf(-dt/TAU_F);
    rDY_[8] =F2_INF-(F2_INF-sf2)*expf(-dt/TAU_F2);

    rDY_[9] = D_INF_new;
    rDY_[10] = R_INF_new;
    rDY_[11] = Xr2_INF_new;


}