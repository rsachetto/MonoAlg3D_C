#include <stdlib.h>

#include "Ohara_Rudy_2011.h"

int celltype = 0;

GET_CELL_MODEL_DATA(init_cell_model_data) {

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

        sv[0] = INITIAL_V; // v
        sv[1] = 7;         // nai
        sv[2] = 7;         // nass
        sv[3] = 145;       // ki
        sv[4] = 145;       // kss
        sv[5] = 1.0e-4;    // cai
        sv[6] = 1.0e-4;    // cass
        sv[7] = 1.2;       // cansr
        sv[8] = 1.2;       // cajsr
        sv[9] = 0;         // m
        sv[10] = 1;        // hf
        sv[11] = 1;        // hs
        sv[12] = 1;        // j
        sv[13] = 1;        // hsp
        sv[14] = 1;        // jp
        sv[15] = 0;        // mL
        sv[16] = 1;        // hL
        sv[17] = 1;        // hLp
        sv[18] = 0;        // a
        sv[19] = 1;        // iF
        sv[20] = 1;        // iS
        sv[21] = 0;        // ap
        sv[22] = 1;        // iFp
        sv[23] = 1;        // iSp
        sv[24] = 0;        // d
        sv[25] = 1;        // ff
        sv[26] = 1;        // fs
        sv[27] = 1;        // fcaf
        sv[28] = 1;        // fcas
        sv[29] = 1;        // jca
        sv[30] = 0;        // nca
        sv[31] = 1;        // ffp
        sv[32] = 1;        // fcafp
        sv[33] = 0;        // xrf
        sv[34] = 0;        // xrs
        sv[35] = 0;        // xs1
        sv[36] = 0;        // xs2
        sv[37] = 1;        // xk1
        sv[38] = 0;        // Jrelnp
        sv[39] = 0;        // Jrelp
        sv[40] = 0;        // CaMKt
    }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;

    int i;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    OMP(parallel for private(sv_id))
    for (i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = (uint32_t )i;

        for (int j = 0; j < num_steps; ++j) {
            solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {


    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt);

    for(int i = 0; i < NEQ; i++)
        sv[i] = rDY[i];
}


void RHS_cpu(const real *sv, real *rDY_, real stim_current,real dt) {

    // State variables
    real v      = sv[0];
    real nai    = sv[1];
    real nass   = sv[2];
    real ki     = sv[3];
    real kss    = sv[4];
    real cai    = sv[5];
    real cass   = sv[6];
    real cansr  = sv[7];
    real cajsr  = sv[8];
    real m      = sv[9];
    real hf     = sv[10];
    real hs     = sv[11];
    real j      = sv[12];
    real hsp    = sv[13];
    real jp     = sv[14];
    real mL     = sv[15];
    real hL     = sv[16];
    real hLp    = sv[17];
    real a      = sv[18];
    real iF     = sv[19];
    real iS     = sv[20];
    real ap     = sv[21];
    real iFp    = sv[22];
    real iSp    = sv[23];
    real d      = sv[24];
    real ff     = sv[25];
    real fs     = sv[26];
    real fcaf   = sv[27];
    real fcas   = sv[28];
    real jca    = sv[29];
    real nca    = sv[30];
    real ffp    = sv[31];
    real fcafp  = sv[32];
    real xrf    = sv[33];
    real xrs    = sv[34];
    real xs1    = sv[35];
    real xs2    = sv[36];
    real xk1    = sv[37];
    real Jrelnp = sv[38];
    real Jrelp  = sv[39];
    real CaMKt  = sv[40];

    //constants
    real const nao=140.0;//extracellular sodium in mM
    real const cao=1.8;//extracellular calcium in mM
    real const ko=5.4;//extracellular potassium in mM

//buffer paramaters
    real const BSRmax=0.047;
    real const KmBSR=0.00087;
    real const BSLmax=1.124;
    real const KmBSL=0.0087;
    real const cmdnmax=0.05;
    real const kmcmdn=0.00238;
    real const trpnmax=0.07;
    real const kmtrpn=0.0005;
    real const csqnmax=10.0f;
    real const kmcsqn=0.8;

//CaMK paramaters
    real const aCaMK=0.05;
    real const bCaMK=0.00068;
    real const CaMKo=0.05;
    real const KmCaM=0.0015;
    real const KmCaMK=0.15;

//physical constants
    real const R=8314.0;
    real const T=310.0;
    real const F=96485.0;

//cell geometry
    real const L=0.01;
    real const rad=0.0011;
    real const vcell=1000.0f*3.14f*rad*rad*L;
    real const Ageo=2.0f*3.14f*rad*rad+2.0f*3.14f*rad*L;
    real const Acap=2*Ageo;
    real const vmyo=0.68f*vcell;
//    real const vmito=0.26*vcell;
//    real const vsr=0.06*vcell;
    real const vnsr=0.0552f*vcell;
    real const vjsr=0.0048f*vcell;
    real const vss=0.02f*vcell;


    real ENa=(R*T/F)*log(nao/nai);
    real EK=(R*T/F)*log(ko/ki);
    real EKs=(R*T/F)*log((ko+0.01833f*nao)/(ki+0.01833f*nai));


    real CaMKb=CaMKo*(1.0f-CaMKt)/(1.0f+KmCaM/cass);
    real CaMKa=CaMKb+CaMKt;
    real vffrt=v*F*F/(R*T);
    real vfrt=v*F/(R*T);

    real mss=1.0f/(1.0f+exp((-(v+39.57f))/9.871f));
    real tm=1.0f/(6.765f*exp((v+11.64f)/34.77f)+8.552f*exp(-(v+77.42f)/5.955f));
    m=mss-(mss-m)*exp(-dt/tm);
    real hss=1.0f/(1+exp((v+82.90f)/6.086f));
    real thf=1.0f/(1.432e-5f*exp(-(v+1.196f)/6.285f)+6.149f*exp((v+0.5096f)/20.27f));
    real ths=1.0f/(0.009794f*exp(-(v+17.95f)/28.05f)+0.3343f*exp((v+5.730f)/56.66f));
    real Ahf=0.99;
    real Ahs=1.0f-Ahf;
    hf=hss-(hss-hf)*exp(-dt/thf);
    hs=hss-(hss-hs)*exp(-dt/ths);
    real h=Ahf*hf+Ahs*hs;
    real jss=hss;
    real tj=2.038f+1.0f/(0.02136f*exp(-(v+100.6f)/8.281f)+0.3052f*exp((v+0.9941f)/38.45f));
    j=jss-(jss-j)*exp(-dt/tj);
    real hssp=1.0f/(1.0f+exp((v+89.1f)/6.086f));
    real thsp=3.0f*ths;
    hsp=hssp-(hssp-hsp)*exp(-dt/thsp);
    real hp=Ahf*hf+Ahs*hsp;
    real tjp=1.46f*tj;
    jp=jss-(jss-jp)*exp(-dt/tjp);
    real GNa=75;
    real fINap=(1.0f/(1.0f+KmCaMK/CaMKa));
    real INa=GNa*(v-ENa)*m*m*m*((1.0f-fINap)*h*j+fINap*hp*jp);

    real mLss=1.0f/(1.0f+exp((-(v+42.85f))/5.264f));
    real tmL=tm;
    mL=mLss-(mLss-mL)*exp(-dt/tmL);
    real hLss=1.0f/(1.0f+exp((v+87.61f)/7.488f));
    real thL=200.0;
    hL=hLss-(hLss-hL)*exp(-dt/thL);
    real hLssp=1.0f/(1.0f+exp((v+93.81f)/7.488f));
    real thLp=3.0f*thL;
    hLp=hLssp-(hLssp-hLp)*exp(-dt/thLp);
    real GNaL=0.0075;
    if (celltype==1)
    {
        GNaL*=0.6;
    }
    real fINaLp=(1.0f/(1.0f+KmCaMK/CaMKa));
    real INaL=GNaL*(v-ENa)*mL*((1.0f-fINaLp)*hL+fINaLp*hLp);

    real ass=1.0f/(1.0f+exp((-(v-14.34f))/14.82f));
    real ta=1.0515f/(1.0f/(1.2089f*(1.0f+exp(-(v-18.4099f)/29.3814f)))+3.5f/(1.0f+exp((v+100.0f)/29.3814f)));
    a=ass-(ass-a)*exp(-dt/ta);
    real iss=1.0f/(1.0f+exp((v+43.94f)/5.711f));
    real delta_epi;
    if (celltype==1)
    {
        delta_epi=1.0f-(0.95f/(1.0f+exp((v+70.0f)/5.0f)));
    }
    else
    {
        delta_epi=1.0f;
    }
    real tiF=4.562f+1.0f/(0.3933f*exp((-(v+100.0f))/100.0f)+0.08004f*exp((v+50.0f)/16.59f));
    real tiS=23.62f+1.0f/(0.001416f*exp((-(v+96.52f))/59.05f)+1.780e-8f*exp((v+114.1f)/8.079f));
    tiF*=delta_epi;
    tiS*=delta_epi;
    real AiF=1.0f/(1.0f+exp((v-213.6f)/151.2f));
    real AiS=1.0f-AiF;
    iF=iss-(iss-iF)*exp(-dt/tiF);
    iS=iss-(iss-iS)*exp(-dt/tiS);
    real i=AiF*iF+AiS*iS;
    real assp=1.0f/(1.0f+exp((-(v-24.34f))/14.82f));
    ap=assp-(assp-ap)*exp(-dt/ta);
    real dti_develop=1.354f+1.0e-4f/(exp((v-167.4f)/15.89f)+exp(-(v-12.23f)/0.2154f));
    real dti_recover=1.0f-0.5f/(1.0f+exp((v+70.0f)/20.0f));
    real tiFp=dti_develop*dti_recover*tiF;
    real tiSp=dti_develop*dti_recover*tiS;
    iFp=iss-(iss-iFp)*exp(-dt/tiFp);
    iSp=iss-(iss-iSp)*exp(-dt/tiSp);
    real ip=AiF*iFp+AiS*iSp;
    real Gto=0.02;
    if (celltype==1)
    {
        Gto*=4.0;
    }
    if (celltype==2)
    {
        Gto*=4.0;
    }
    real fItop=(1.0f/(1.0f+KmCaMK/CaMKa));
    real Ito = Gto*(v-EK)*((1.0f-fItop)*a*i+fItop*ap*ip);

    real dss=1.0f/(1.0f+exp((-(v+3.940f))/4.230f));
    real td=0.6f+1.0f/(exp(-0.05f*(v+6.0f))+exp(0.09f*(v+14.0f)));
    d=dss-(dss-d)*exp(-dt/td);
    real fss=1.0f/(1.0f+exp((v+19.58f)/3.696f));
    real tff=7.0f+1.0f/(0.0045f*exp(-(v+20.0f)/10.0f)+0.0045f*exp((v+20.0f)/10.0f));
    real tfs=1000.0f+1.0f/(0.000035f*exp(-(v+5.0f)/4.0f)+0.000035f*exp((v+5.0f)/6.0f));
    real Aff=0.6;
    real Afs=1.0f-Aff;
    ff=fss-(fss-ff)*exp(-dt/tff);
    fs=fss-(fss-fs)*exp(-dt/tfs);
    real f=Aff*ff+Afs*fs;
    real fcass=fss;
    real tfcaf=7.0f+1.0f/(0.04f*exp(-(v-4.0f)/7.0f)+0.04f*exp((v-4.0f)/7.0f));
    real tfcas=100.0f+1.0f/(0.00012f*exp(-v/3.0f)+0.00012f*exp(v/7.0f));
    real Afcaf=0.3f+0.6f/(1.0f+exp((v-10.0f)/10.0f));
    real Afcas=1.0f-Afcaf;
    fcaf=fcass-(fcass-fcaf)*exp(-dt/tfcaf);
    fcas=fcass-(fcass-fcas)*exp(-dt/tfcas);
    real fca=Afcaf*fcaf+Afcas*fcas;
    real tjca=75.0;
    jca=fcass-(fcass-jca)*exp(-dt/tjca);
    real tffp=2.5f*tff;
    ffp=fss-(fss-ffp)*exp(-dt/tffp);
    real fp=Aff*ffp+Afs*fs;
    real tfcafp=2.5f*tfcaf;
    fcafp=fcass-(fcass-fcafp)*exp(-dt/tfcafp);
    real fcap=Afcaf*fcafp+Afcas*fcas;
    real Kmn=0.002;
    real k2n=1000.0;
    real km2n=jca*1.0f;
    real anca=1.0f/(k2n/km2n+pow(1.0f+Kmn/cass,4.0));
    nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
    real PhiCaL=4.0f*vffrt*(cass*exp(2.0f*vfrt)-0.341f*cao)/(exp(2.0f*vfrt)-1.0f);
    real PhiCaNa=1.0f*vffrt*(0.75f*nass*exp(1.0f*vfrt)-0.75f*nao)/(exp(1.0f*vfrt)-1.0f);
    real PhiCaK=1.0f*vffrt*(0.75f*kss*exp(1.0f*vfrt)-0.75f*ko)/(exp(1.0f*vfrt)-1.0f);
    real zca=2.0f;
    real PCa=0.0001;
    if (celltype==1)
    {
        PCa*=1.2;
    }
    if (celltype==2)
    {
        PCa*=2.5;
    }
    real PCap=1.1f*PCa;
    real PCaNa=0.00125f*PCa;
    real PCaK=3.574e-4f*PCa;
    real PCaNap=0.00125f*PCap;
    real PCaKp=3.574e-4f*PCap;
    real fICaLp=(1.0f/(1.0f+KmCaMK/CaMKa));
    real ICaL=(1.0f-fICaLp)*PCa*PhiCaL*d*(f*(1.0f-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0f-nca)+jca*fcap*nca);
    real ICaNa=(1.0f-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0f-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0f-nca)+jca*fcap*nca);
    real ICaK=(1.0f-fICaLp)*PCaK*PhiCaK*d*(f*(1.0f-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0f-nca)+jca*fcap*nca);

    real xrss=1.0f/(1.0f+exp((-(v+8.337f))/6.789f));
    real txrf=12.98f+1.0f/(0.3652f*exp((v-31.66f)/3.869f)+4.123e-5f*exp((-(v-47.78f))/20.38f));
    real txrs=1.865f+1.0f/(0.06629f*exp((v-34.70f)/7.355f)+1.128e-5f*exp((-(v-29.74f))/25.94f));
    real Axrf=1.0f/(1.0f+exp((v+54.81f)/38.21f));
    real Axrs=1.0f-Axrf;
    xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
    xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
    real xr=Axrf*xrf+Axrs*xrs;
    real rkr=1.0f/(1.0f+exp((v+55.0f)/75.0f))*1.0f/(1.0f+exp((v-10.0f)/30.0f));
    real GKr=0.046;
    if (celltype==1)
    {
        GKr*=1.3;
    }
    if (celltype==2)
    {
        GKr*=0.8;
    }
    real IKr = GKr*sqrt(ko/5.4f)*xr*rkr*(v-EK);

    real xs1ss=1.0f/(1.0f+exp((-(v+11.60f))/8.932f));
    real txs1=817.3f+1.0f/(2.326e-4f*exp((v+48.28f)/17.80f)+0.001292f*exp((-(v+210.0f))/230.0f));
    xs1=xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
    real xs2ss=xs1ss;
    real txs2=1.0f/(0.01f*exp((v-50.0f)/20.0f)+0.0193f*exp((-(v+66.54f))/31.0f));
    xs2=xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
    real KsCa=1.0f+0.6f/(1.0f+pow(3.8e-5f/cai,1.4f));
    real GKs=0.0034;
    if (celltype==1)
    {
        GKs*=1.4;
    }
    real  IKs=GKs*KsCa*xs1*xs2*(v-EKs);

    real xk1ss=1.0f/(1.0f+exp(-(v+2.5538f*ko+144.59f)/(1.5692f*ko+3.8115f)));
    real txk1=122.2f/(exp((-(v+127.2f))/20.36f)+exp((v+236.8f)/69.33f));
    xk1=xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
    real rk1=1.0f/(1.0f+exp((v+105.8f-2.6f*ko)/9.493f));
    real GK1=0.1908;
    if (celltype==1)
    {
        GK1*=1.2;
    }
    if (celltype==2)
    {
        GK1*=1.3;
    }
    real  IK1=GK1*sqrt(ko)*rk1*xk1*(v-EK);

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
    real h3=1.0f/h1;
    real h4=1.0f+nai/kna1*(1+nai/kna2);
    real h5=nai*nai/(h4*kna1*kna2);
    real h6=1.0f/h4;
    real h7=1.0f+nao/kna3*(1.0f+1.0f/hna);
    real h8=nao/(kna3*hna*h7);
    real h9=1.0f/h7;
    real h10=kasymm+1.0f+nao/kna1*(1.0f+nao/kna2);
    real h11=nao*nao/(h10*kna1*kna2);
    real h12=1.0f/h10;
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
    real allo=1.0f/(1.0f+pow(KmCaAct/cai,2.0f));
    real zna=1.0f;
    real JncxNa=3.0f*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    real JncxCa=E2*k2-E1*k1;
    real Gncx=0.0008;
    if (celltype==1)
    {
        Gncx*=1.1;
    }
    if (celltype==2)
    {
        Gncx*=1.4;
    }

    real INaCa_i=0.8f*Gncx*allo*(zna*JncxNa+zca*JncxCa);

    h1=1+nass/kna3*(1+hna);
    h2=(nass*hna)/(kna3*h1);
    h3=1.0f/h1;
    h4=1.0f+nass/kna1*(1+nass/kna2);
    h5=nass*nass/(h4*kna1*kna2);
    h6=1.0f/h4;
    h7=1.0f+nao/kna3*(1.0f+1.0f/hna);
    h8=nao/(kna3*hna*h7);
    h9=1.0f/h7;
    h10=kasymm+1.0f+nao/kna1*(1+nao/kna2);
    h11=nao*nao/(h10*kna1*kna2);
    h12=1.0f/h10;
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
    KmCaAct=150.0e-6;
    allo=1.0f/(1.0f+pow(KmCaAct/cass,2.0f));
    JncxNa=3.0f*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    JncxCa=E2*k2-E1*k1;
    real INaCa_ss=0.2f*Gncx*allo*(zna*JncxNa+zca*JncxCa);

    real INaCa=INaCa_i+INaCa_ss;

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
    real delta=-0.1550f;
    real Knai=Knai0*exp((delta*v*F)/(3.0f*R*T));
    real Knao=Knao0*exp(((1.0f-delta)*v*F)/(3.0f*R*T));
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
    real P=eP/(1.0f+H/Khp+nai/Knap+ki/Kxkur);
    real a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0f+nai/Knai,3.0)+pow(1.0f+ki/Kki,2.0f)-1.0f);
    real b1=k1m*MgADP;
    real a2=k2p;
    real b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0f+nao/Knao,3.0)+pow(1.0f+ko/Kko,2.0f)-1.0f);
    real a3=(k3p*pow(ko/Kko,2.0f))/(pow(1.0f+nao/Knao,3.0)+pow(1.0f+ko/Kko,2.0f)-1.0f);
    real b3=(k3m*P*H)/(1.0f+MgATP/Kmgatp);
    real a4=(k4p*MgATP/Kmgatp)/(1.0f+MgATP/Kmgatp);
    real b4=(k4m*pow(ki/Kki,2.0f))/(pow(1.0f+nai/Knai,3.0)+pow(1.0f+ki/Kki,2.0f)-1.0f);
    x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
    x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
    x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
    x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
    E1=x1/(x1+x2+x3+x4);
    E2=x2/(x1+x2+x3+x4);
    E3=x3/(x1+x2+x3+x4);
    E4=x4/(x1+x2+x3+x4);
    real zk=1.0f;
    real JnakNa=3.0f*(E1*a3-E2*b3);
    real JnakK=2.0f*(E4*b1-E3*a1);
    real Pnak=30;
    if (celltype==1)
    {
        Pnak*=0.9;
    }
    if (celltype==2)
    {
        Pnak*=0.7;
    }
    real INaK=Pnak*(zna*JnakNa+zk*JnakK);

    real xkb=1.0f/(1.0f+exp(-(v-14.48f)/18.34f));
    real GKb=0.003;
    if (celltype==1)
    {
        GKb*=0.6;
    }
    real IKb=GKb*xkb*(v-EK);

    real PNab=3.75e-10;
    real INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0f);

    real PCab=2.5e-8;
    real ICab=PCab*4.0f*vffrt*(cai*exp(2.0f*vfrt)-0.341f*cao)/(exp(2.0f*vfrt)-1.0f);

    real GpCa=0.0005;
    real IpCa=GpCa*cai/(0.0005f+cai);

    //voltage()
    v+=-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+stim_current);

    CaMKb = CaMKo*(1.0f-CaMKt)/(1.0f+KmCaM/cass);
    CaMKa = CaMKb+CaMKt;
    CaMKt+=dt*(aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt);

    real JdiffNa=(nass-nai)/2.0f;
    real JdiffK=(kss-ki)/2.0f;
    real Jdiff=(cass-cai)/0.2f;

    real bt=4.75;
    real a_rel=0.5f*bt;
    real Jrel_inf=a_rel*(-ICaL)/(1.0f+pow(1.5f/cajsr,8.0));
    if (celltype==2)
    {
        Jrel_inf*=1.7;
    }
    real tau_rel=bt/(1.0f+0.0123f/cajsr);
    if (tau_rel<0.005)
    {
        tau_rel=0.005;
    }
    Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
    real btp=1.25f*bt;
    real a_relp=0.5f*btp;
    real Jrel_infp=a_relp*(-ICaL)/(1.0f+pow(1.5f/cajsr,8.0f));
    if (celltype==2)
    {
        Jrel_infp*=1.7;
    }
    real tau_relp=btp/(1.0f+0.0123f/cajsr);
    if (tau_relp<0.005f)
    {
        tau_relp=0.005f;
    }
    Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp);
    real fJrelp=(1.0f/(1.0f+KmCaMK/CaMKa));
    real Jrel=(1.0f-fJrelp)*Jrelnp+fJrelp*Jrelp;

    real Jupnp=0.004375f*cai/(cai+0.00092f);
    real Jupp=2.75f*0.004375f*cai/(cai+0.00092f-0.00017f);
    if (celltype==1)
    {
        Jupnp*=1.3;
        Jupp*=1.3;
    }
    real fJupp=(1.0f/(1.0f+KmCaMK/CaMKa));
    real Jleak=0.0039375f*cansr/15.0f;
    real Jup=(1.0f-fJupp)*Jupnp+fJupp*Jupp-Jleak;

    real Jtr=(cansr-cajsr)/100.0f;

    nai+=dt*(-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);
    nass+=dt*(-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);

    ki+=dt*(-(Ito+IKr+IKs+IK1+IKb+stim_current-2.0f*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);
    kss+=dt*(-(ICaK)*Acap/(F*vss)-JdiffK);

    real Bcai;
    if (celltype==1)
    {
        Bcai=1.0f/(1.0f+1.3f*cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0f)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0f));
    }
    else
    {
        Bcai=1.0f/(1.0f+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0f)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0f));
    }
    cai+=dt*(Bcai*(-(IpCa+ICab-2.0f*INaCa_i)*Acap/(2.0f*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));

    real Bcass=1.0f/(1.0f+BSRmax*KmBSR/pow(KmBSR+cass,2.0f)+BSLmax*KmBSL/pow(KmBSL+cass,2.0f));
    cass+=dt*(Bcass*(-(ICaL-2.0f*INaCa_ss)*Acap/(2.0f*F*vss)+Jrel*vjsr/vss-Jdiff));

    cansr+=dt*(Jup-Jtr*vjsr/vnsr);

    real Bcajsr=1.0f/(1.0f+csqnmax*kmcsqn/pow(kmcsqn+cajsr,2.0f));
    cajsr+=dt*(Bcajsr*(Jtr-Jrel));


    rDY_[0]  =  v;
    rDY_[1]  =  nai;
    rDY_[2]  =  nass;
    rDY_[3]  =  ki;
    rDY_[4]  =  kss;
    rDY_[5]  =  cai;
    rDY_[6]  =  cass;
    rDY_[7]  =  cansr;
    rDY_[8]  =  cajsr;
    rDY_[9]  =  m;
    rDY_[10] =  hf;
    rDY_[11] =  hs;
    rDY_[12] =  j;
    rDY_[13] =  hsp;
    rDY_[14] =  jp;
    rDY_[15] =  mL;
    rDY_[16] =  hL;
    rDY_[17] =  hLp;
    rDY_[18] =  a;
    rDY_[19] =  iF;
    rDY_[20] =  iS;
    rDY_[21] =  ap;
    rDY_[22] =  iFp;
    rDY_[23] =  iSp;
    rDY_[24] =  d;
    rDY_[25] =  ff;
    rDY_[26] =  fs;
    rDY_[27] =  fcaf;
    rDY_[28] =  fcas;
    rDY_[29] =  jca;
    rDY_[30] =  nca;
    rDY_[31] =  ffp;
    rDY_[32] =  fcafp;
    rDY_[33] =  xrf;
    rDY_[34] =  xrs;
    rDY_[35] =  xs1;
    rDY_[36] =  xs2;
    rDY_[37] =  xk1;
    rDY_[38] =  Jrelnp;
    rDY_[39] =  Jrelp;
    rDY_[40] =  CaMKt;


}
