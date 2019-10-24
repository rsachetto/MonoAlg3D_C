
#include "model_common.h"
#include <assert.h>
#include <stdlib.h>
#include "ten_tusscher_3_sensibility.h"

#define ENDO

GET_CELL_MODEL_DATA(init_cell_model_data) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;

}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

	// Default initial conditions
    if (extra_data == NULL)
    {
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
    else
    {
        real *initial_conditions = ((real*)extra_data) + 7; //pointer
        real *fibrosis = ((real*)extra_data) + 7 + 12;      //pointer

        sv[0] = initial_conditions[0];   // V;       millivolt
		sv[1] = initial_conditions[1]; //M
		sv[2] = initial_conditions[2]; //H
		sv[3] = initial_conditions[3]; //J
		sv[4] = initial_conditions[4]; //Xr1
		sv[5] = initial_conditions[5]; //Xs
		sv[6] = initial_conditions[6]; //S
		sv[7] = initial_conditions[7]; //F
		sv[8] = initial_conditions[8]; //F2
		sv[9] = initial_conditions[9]; //D_INF
		sv[10] = initial_conditions[10]; //R_INF
		sv[11] = initial_conditions[11]; //Xr2_INF}
    }
}

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    uint32_t sv_id;
    real *fibrosis;

    // Default values for a healthy cell ///////////
    real atpi = 6.8f;
    real Ko = 5.4f;
    real Ki = 138.3f;
    real Vm_change = 0.0;
    real GNa_multiplicator = 1.0f;
    real GCa_multiplicator = 1.0f;
    real INaCa_multiplicator = 1.0f;
    
    real V_0 = -86.2f;   // V;       millivolt
	real M_0 = 0.0f; //M
	real H_0 = 0.75; //H
	real J_0 = 0.75; //J
	real Xr1_0 = 0.0f; //Xr1
	real Xs_0 = 0.0f; //Xs
	real S_0 = 1.0f; //S
	real F_0 = 1.0f; //F
	real F2_0 = 1.0f; //F2
	real D_inf_0 = 0.0; //D_INF
	real R_inf_0 = 0.0; //R_INF
	real Xr2_inf_0 = 0.0; //Xr2_INF}
    ////////////////////////////////////

    int num_extra_parameters = 7;
    int num_initial_conditions = 12;
    size_t extra_parameters_size = (num_extra_parameters+num_initial_conditions)*sizeof(real);

    if(extra_data) 
    {
        fibrosis = ((real*)extra_data) + num_extra_parameters + num_initial_conditions; //pointer
    }
    else 
    {
        extra_data = malloc(extra_parameters_size);
        ((real*)extra_data)[0] = atpi;
        ((real*)extra_data)[1] = Ko;
        ((real*)extra_data)[2] = Ki;
        ((real*)extra_data)[3] = Vm_change;
        ((real*)extra_data)[4] = GNa_multiplicator;
        ((real*)extra_data)[5] = GCa_multiplicator;
        ((real*)extra_data)[6] = INaCa_multiplicator;
        
        ((real*)extra_data)[7] = V_0;
        ((real*)extra_data)[8] = M_0;
        ((real*)extra_data)[9] = H_0;
        ((real*)extra_data)[10] = J_0;
        ((real*)extra_data)[11] = Xr1_0;
        ((real*)extra_data)[12] = Xs_0;
        ((real*)extra_data)[13] = S_0;
        ((real*)extra_data)[14] = F_0;
        ((real*)extra_data)[15] = F2_0;
        ((real*)extra_data)[16] = D_inf_0;
        ((real*)extra_data)[17] = R_inf_0;
        ((real*)extra_data)[18] = Xr2_inf_0;

        fibrosis = calloc(num_cells_to_solve, sizeof(real));
    }

    int i;

	#pragma omp parallel for private(sv_id)
    for (i = 0; i < num_cells_to_solve; i++) 
    {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        for (int j = 0; j < num_steps; ++j) 
        {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], fibrosis[i], extra_data);
        }
    }

    if(extra_data == NULL) free(fibrosis);
}


void solve_model_ode_cpu(real dt, real *sv, real stim_current, real fibrosis, real *extra_parameters)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, fibrosis, extra_parameters);

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


void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real fibrosis, real *extra_parameters) {

    //fibrosis = 0 means that the cell is fibrotic, 1 is not fibrotic. Anything between 0 and 1 means border zone
        
    //THIS IS THE STATE VECTOR THAT WE NEED TO SAVE IN THE STEADY STATE
    const real svolt    = sv[0];
    const real sm       = sv[1];
    const real sh       = sv[2];
    const real sj       = sv[3];
    const real sxr1     = sv[4];
    const real sxs      = sv[5];
    const real ss       = sv[6];
    const real sf       = sv[7];
    const real sf2      = sv[8];
    const real D_INF    = sv[9];
    const real R_INF    = sv[10];
    const real Xr2_INF  = sv[11];

    const real natp = 0.24;          // K dependence of ATP-sensitive K current
    const real nicholsarea = 0.00005; // Nichol's areas (cm^2)
    const real hatp = 2;             // Hill coefficient

    //Linear changing of atpi depending on the fibrosis and distance from the center of the scar (only for border zone cells)
    real atpi = extra_parameters[0];
    real atpi_change = 6.8f - atpi;
    atpi = atpi + atpi_change*fibrosis;

    //Extracellular potassium concentration was elevated
    //from its default value of 5.4 mM to values between 6.0 and 8.0 mM
    //Ref: A Comparison of Two Models of Human Ventricular Tissue: Simulated Ischemia and Re-entry
    real Ko = extra_parameters[1];
    real Ko_change  = 5.4f - Ko;
    Ko = Ko + Ko_change*fibrosis;

    real Ki = extra_parameters[2];
    real Ki_change  = 138.3 - Ki;
    Ki = Ki + Ki_change*fibrosis;  

    real Vm_modifier = extra_parameters[3];
    Vm_modifier = Vm_modifier - Vm_modifier*fibrosis;

    real GNa_multplicator = extra_parameters[4];
    real GNa_multplicator_change  = 1.0f - GNa_multplicator;
    GNa_multplicator = GNa_multplicator + GNa_multplicator_change*fibrosis;

    real GCaL_multplicator = extra_parameters[5];
    real GCaL_multplicator_change  = 1.0f - GCaL_multplicator;
    GCaL_multplicator = GCaL_multplicator + GCaL_multplicator_change*fibrosis;

    real INaCa_multplicator = extra_parameters[6];
    real INaCa_multplicator_change  = 1.0f - INaCa_multplicator;
    INaCa_multplicator = INaCa_multplicator + INaCa_multplicator_change*fibrosis;
     

    //real katp = 0.306;
    //Ref: A Comparison of Two Models of Human Ventricular Tissue: Simulated Ischaemia and Re-entry
    //real katp = 0.306;
    const real katp = -0.0942857142857*atpi + 0.683142857143; //Ref: A Comparison of Two Models of Human Ventricular Tissue: Simulated Ischaemia and Re-entry

    const real patp =  1/(1 + pow((atpi/katp),hatp));
    const real gkatp    =  0.000195/nicholsarea;
    const real gkbaratp =  gkatp*patp*pow((Ko/5.4),natp);

    const real katp2= 1.4;
    const real hatp2 = 2.6;
    const real pcal = 1.0/(1.0 + pow((katp2/atpi),hatp2));


    const real Cao=2.0;
    const real Nao=140.0;
    const real Cai=0.00007;
    const real Nai=7.67;

//Constants
    const real R=8314.472;
    const real F=96485.3415;
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
    const real GNa=14.838*GNa_multplicator; //ACIDOSIS
//Parameters for IbNa
    const real GbNa=0.00029;
//Parameters for INaK
    const real KmK=1.0;
    const real KmNa=40.0;
    const real knak=2.724;
//Parameters for ICaL
    const real GCaL=0.2786*pcal*GCaL_multplicator; //ACIDOSIS
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

    const real Ek=RTONF*(log((Ko/Ki)));
    const real Ena=RTONF*(log((Nao/Nai)));
    const real Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    const real Eca=0.5*RTONF*(log((Cao/Cai)));
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
    Ak1=0.1/(1.+exp(0.06*(svolt-Ek-200)));
    Bk1=(3.*exp(0.0002*(svolt-Ek+100))+
         exp(0.1*(svolt-Ek-10)))/(1.+exp(-0.5*(svolt-Ek)));
    rec_iK1=Ak1/(Ak1+Bk1);
    rec_iNaK=(1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
    rec_ipK=1./(1.+exp((25-svolt)/5.98));


    //Compute currents
    INa=GNa*sm*sm*sm*sh*sj*((svolt-Vm_modifier)-Ena); //ACIDOSIS
    ICaL=GCaL*D_INF*sf*sf2*((svolt-Vm_modifier)-60); //ACIDOSIS
    Ito=Gto*R_INF*ss*(svolt-Ek);
    IKr=Gkr*sqrt(Ko/5.4)*sxr1*Xr2_INF*(svolt-Ek);
    IKs=Gks*sxs*sxs*(svolt-Eks);
    IK1=GK1*rec_iK1*(svolt-Ek);
    INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
          (1./(1+ksat*exp((n-1)*svolt*F/(R*T))))*
          (exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
           exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);

    INaCa = INaCa*INaCa_multplicator; //ACIDOSIS

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
