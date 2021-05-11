#include "../../gpu_utils/gpu_utils.h"
#include <stddef.h>
#include <stdint.h>

#include "Ohara_Rudy_2011.h"
__device__ int celltype_gpu = 0;

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    log_info("Using Ohara Rudy 2011 GPU model\n");

    uint32_t num_volumes = solver->original_num_cells;

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(solver->sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}


extern "C" SOLVE_MODEL_ODES(solve_model_odes_gpu) {


    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    //the array cells to solve is passed when we are using and adaptive mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }
    solve_gpu <<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_inital_conditions (real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        *((real * )((char *) sv + pitch * 0) + threadID) =   INITIAL_V;      //v
        *((real * )((char *) sv + pitch * 1) + threadID) =   7;              //nai
        *((real * )((char *) sv + pitch * 2) + threadID) =   7;              //nass
        *((real * )((char *) sv + pitch * 3) + threadID) =   145;            //ki
        *((real * )((char *) sv + pitch * 4) + threadID) =   145;            //kss
        *((real * )((char *) sv + pitch * 5) + threadID) =   1.0e-4;         //cai
        *((real * )((char *) sv + pitch * 6) + threadID) =   1.0e-4;         //cass
        *((real * )((char *) sv + pitch * 7) + threadID) =   1.2;            //cansr
        *((real * )((char *) sv + pitch * 8) + threadID) =   1.2;            //cajsr
        *((real * )((char *) sv + pitch * 9) + threadID) =   0;              //m
        *((real * )((char *) sv + pitch * 10) + threadID) =  1;              //hf
        *((real * )((char *) sv + pitch * 11) + threadID) =  1;              //hs
        *((real * )((char *) sv + pitch * 12) + threadID) =  1;              //j
        *((real * )((char *) sv + pitch * 13) + threadID) =  1;              //hsp
        *((real * )((char *) sv + pitch * 14) + threadID) =  1;              //jp
        *((real * )((char *) sv + pitch * 15) + threadID) =  0;              //mL
        *((real * )((char *) sv + pitch * 16) + threadID) =  1;              //hL
        *((real * )((char *) sv + pitch * 17) + threadID) =  1;              //hLp
        *((real * )((char *) sv + pitch * 18) + threadID) =  0;              //a
        *((real * )((char *) sv + pitch * 19) + threadID) =  1;              //iF
        *((real * )((char *) sv + pitch * 20) + threadID) =  1;              //iS
        *((real * )((char *) sv + pitch * 21) + threadID) =  0;              //ap
        *((real * )((char *) sv + pitch * 22) + threadID) =  1;              //iFp
        *((real * )((char *) sv + pitch * 23) + threadID) =  1;              //iSp
        *((real * )((char *) sv + pitch * 24) + threadID) =  0;              //d
        *((real * )((char *) sv + pitch * 25) + threadID) =  1;              //ff
        *((real * )((char *) sv + pitch * 26) + threadID) =  1;              //fs
        *((real * )((char *) sv + pitch * 27) + threadID) =  1;              //fcaf
        *((real * )((char *) sv + pitch * 28) + threadID) =  1;              //fcas
        *((real * )((char *) sv + pitch * 29) + threadID) =  1;              //jca
        *((real * )((char *) sv + pitch * 30) + threadID) =  0;              //nca
        *((real * )((char *) sv + pitch * 31) + threadID) =  1;              //ffp
        *((real * )((char *) sv + pitch * 32) + threadID) =  1;              //fcafp
        *((real * )((char *) sv + pitch * 33) + threadID) =  0;              //xrf
        *((real * )((char *) sv + pitch * 34) + threadID) =  0;              //xrs
        *((real * )((char *) sv + pitch * 35) + threadID) =  0;              //xs1
        *((real * )((char *) sv + pitch * 36) + threadID) =  0;              //xs2
        *((real * )((char *) sv + pitch * 37) + threadID) =  1;              //xk1
        *((real * )((char *) sv + pitch * 38) + threadID) =  0;              //Jrelnp
        *((real * )((char *) sv + pitch * 39) + threadID) =  0;              //Jrelp
        *((real * )((char *) sv + pitch * 40) + threadID) =  0;              //CaMKt

    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu (real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve,
                           uint32_t num_cells_to_solve, int num_steps) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if (threadID < num_cells_to_solve) {
        if (cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            RHS_gpu (sv, rDY, stim_currents[threadID], sv_id, dt);

            for (int i = 0; i < NEQ; i++) {
                *((real *)((char *)sv + pitch * i) + sv_id) = rDY[i];
            }
        }
    }
}

inline __device__ void RHS_gpu (real *sv_, real *rDY_, real stim_current, int threadID_, real dt) {

    // State variables
    real v      = *((real *)((char *)sv_ + pitch * 0) + threadID_);
    real nai    = *((real *)((char *)sv_ + pitch * 1) + threadID_);
    real nass   = *((real *)((char *)sv_ + pitch * 2) + threadID_);
    real ki     = *((real *)((char *)sv_ + pitch * 3) + threadID_);
    real kss    = *((real *)((char *)sv_ + pitch * 4) + threadID_);
    real cai    = *((real *)((char *)sv_ + pitch * 5) + threadID_);
    real cass   = *((real *)((char *)sv_ + pitch * 6) + threadID_);
    real cansr  = *((real *)((char *)sv_ + pitch * 7) + threadID_);
    real cajsr  = *((real *)((char *)sv_ + pitch * 8) + threadID_);
    real m      = *((real *)((char *)sv_ + pitch * 9) + threadID_);
    real hf     = *((real *)((char *)sv_ + pitch * 10) + threadID_);
    real hs     = *((real *)((char *)sv_ + pitch * 11) + threadID_);
    real j      = *((real *)((char *)sv_ + pitch * 12) + threadID_);
    real hsp    = *((real *)((char *)sv_ + pitch * 13) + threadID_);
    real jp     = *((real *)((char *)sv_ + pitch * 14) + threadID_);
    real mL     = *((real *)((char *)sv_ + pitch * 15) + threadID_);
    real hL     = *((real *)((char *)sv_ + pitch * 16) + threadID_);
    real hLp    = *((real *)((char *)sv_ + pitch * 17) + threadID_);
    real a      = *((real *)((char *)sv_ + pitch * 18) + threadID_);
    real iF     = *((real *)((char *)sv_ + pitch * 19) + threadID_);
    real iS     = *((real *)((char *)sv_ + pitch * 20) + threadID_);
    real ap     = *((real *)((char *)sv_ + pitch * 21) + threadID_);
    real iFp    = *((real *)((char *)sv_ + pitch * 22) + threadID_);
    real iSp    = *((real *)((char *)sv_ + pitch * 23) + threadID_);
    real d      = *((real *)((char *)sv_ + pitch * 24) + threadID_);
    real ff     = *((real *)((char *)sv_ + pitch * 25) + threadID_);
    real fs     = *((real *)((char *)sv_ + pitch * 26) + threadID_);
    real fcaf   = *((real *)((char *)sv_ + pitch * 27) + threadID_);
    real fcas   = *((real *)((char *)sv_ + pitch * 28) + threadID_);
    real jca    = *((real *)((char *)sv_ + pitch * 29) + threadID_);
    real nca    = *((real *)((char *)sv_ + pitch * 30) + threadID_);
    real ffp    = *((real *)((char *)sv_ + pitch * 31) + threadID_);
    real fcafp  = *((real *)((char *)sv_ + pitch * 32) + threadID_);
    real xrf    = *((real *)((char *)sv_ + pitch * 33) + threadID_);
    real xrs    = *((real *)((char *)sv_ + pitch * 34) + threadID_);
    real xs1    = *((real *)((char *)sv_ + pitch * 35) + threadID_);
    real xs2    = *((real *)((char *)sv_ + pitch * 36) + threadID_);
    real xk1    = *((real *)((char *)sv_ + pitch * 37) + threadID_);
    real Jrelnp = *((real *)((char *)sv_ + pitch * 38) + threadID_);
    real Jrelp  = *((real *)((char *)sv_ + pitch * 39) + threadID_);
    real CaMKt  = *((real *)((char *)sv_ + pitch * 40) + threadID_);

    //constants
    real const nao=140.0f;//extracellular sodium in mM
    real const cao=1.8f;//extracellular calcium in mM
    real const ko=5.4f;//extracellular potassium in mM

//buffer paramaters
    real const BSRmax=0.047f;
    real const KmBSR=0.00087f;
    real const BSLmax=1.124f;
    real const KmBSL=0.0087f;
    real const cmdnmax=0.05f;
    real const kmcmdn=0.00238f;
    real const trpnmax=0.07f;
    real const kmtrpn=0.0005f;
    real const csqnmax=10.0f;
    real const kmcsqn=0.8f;

//CaMK paramaters
    real const aCaMK=0.05f;
    real const bCaMK=0.00068f;
    real const CaMKo=0.05f;
    real const KmCaM=0.0015f;
    real const KmCaMK=0.15f;

//physical constants
    real const R=8314.0f;
    real const T=310.0f;
    real const F=96485.0f;

//cell geometry
    real const L=0.01f;
    real const rad=0.0011f;
    real const vcell=1000.0f*3.14f*rad*rad*L;
    real const Ageo=2.0f*3.14f*rad*rad+2.0f*3.14f*rad*L;
    real const Acap=2.0f*Ageo;
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
    real thf=1.0f/(1.432e-5*exp(-(v+1.196f)/6.285f)+6.149f*exp((v+0.5096f)/20.27f));
    real ths=1.0f/(0.009794f*exp(-(v+17.95f)/28.05f)+0.3343f*exp((v+5.730f)/56.66f));
    real Ahf=0.99f;
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
    real GNa=75.0f;
    real fINap=(1.0f/(1.0f+KmCaMK/CaMKa));
    real INa=GNa*(v-ENa)*m*m*m*((1.0f-fINap)*h*j+fINap*hp*jp);

    real mLss=1.0/(1.0+exp((-(v+42.85f))/5.264f));
    real tmL=tm;
    mL=mLss-(mLss-mL)*exp(-dt/tmL);
    real hLss=1.0f/(1.0f+exp((v+87.61f)/7.488f));
    real thL=200.0f;
    hL=hLss-(hLss-hL)*exp(-dt/thL);
    real hLssp=1.0f/(1.0f+exp((v+93.81f)/7.488f));
    real thLp=3.0f*thL;
    hLp=hLssp-(hLssp-hLp)*exp(-dt/thLp);
    real GNaL=0.0075f;
    if (celltype_gpu==1)
    {
        GNaL*=0.6f;
    }
    real fINaLp=(1.0f/(1.0f+KmCaMK/CaMKa));
    real INaL=GNaL*(v-ENa)*mL*((1.0f-fINaLp)*hL+fINaLp*hLp);

    real ass=1.0f/(1.0f+exp((-(v-14.34f))/14.82f));
    real ta=1.0515f/(1.0f/(1.2089f*(1.0f+exp(-(v-18.4099f)/29.3814f)))+3.5f/(1.0f+exp((v+100.0f)/29.3814f)));
    a=ass-(ass-a)*exp(-dt/ta);
    real iss=1.0f/(1.0f+exp((v+43.94f)/5.711f));
    real delta_epi;
    if (celltype_gpu==1)
    {
        delta_epi=1.0f-(0.95f/(1.0f+exp((v+70.0f)/5.0f)));
    }
    else
    {
        delta_epi=1.0f;
    }
    real tiF=4.562f+1.0f/(0.3933f*exp((-(v+100.0f))/100.0f)+0.08004f*exp((v+50.0f)/16.59f));
    real tiS=23.62f+1.0f/(0.001416f*exp((-(v+96.52f))/59.05f)+1.780e-8*exp((v+114.1f)/8.079f));
    tiF*=delta_epi;
    tiS*=delta_epi;
    real AiF=1.0f/(1.0f+exp((v-213.6f)/151.2f));
    real AiS=1.0f-AiF;
    iF=iss-(iss-iF)*exp(-dt/tiF);
    iS=iss-(iss-iS)*exp(-dt/tiS);
    real i=AiF*iF+AiS*iS;
    real assp=1.0f/(1.0f+exp((-(v-24.34f))/14.82f));
    ap=assp-(assp-ap)*exp(-dt/ta);
    real dti_develop=1.354f+1.0e-4/(exp((v-167.4f)/15.89f)+exp(-(v-12.23f)/0.2154f));
    real dti_recover=1.0f-0.5/(1.0f+exp((v+70.0f)/20.0f));
    real tiFp=dti_develop*dti_recover*tiF;
    real tiSp=dti_develop*dti_recover*tiS;
    iFp=iss-(iss-iFp)*exp(-dt/tiFp);
    iSp=iss-(iss-iSp)*exp(-dt/tiSp);
    real ip=AiF*iFp+AiS*iSp;
    real Gto=0.02f;
    if (celltype_gpu==1)
    {
        Gto*=4.0f;
    }
    if (celltype_gpu==2)
    {
        Gto*=4.0f;
    }
    real fItop=(1.0f/(1.0f+KmCaMK/CaMKa));
    real Ito = Gto*(v-EK)*((1.0f-fItop)*a*i+fItop*ap*ip);

    real dss=1.0f/(1.0f+exp((-(v+3.940f))/4.230f));
    real td=0.6f+1.0f/(exp(-0.05f*(v+6.0f))+exp(0.09f*(v+14.0f)));
    d=dss-(dss-d)*exp(-dt/td);
    real fss=1.0f/(1.0f+exp((v+19.58)/3.696));
    real tff=7.0f+1.0f/(0.0045f*exp(-(v+20.0f)/10.0f)+0.0045f*exp((v+20.0f)/10.0f));
    real tfs=1000.0f+1.0f/(0.000035f*exp(-(v+5.0f)/4.0f)+0.000035f*exp((v+5.0f)/6.0));
    real Aff=0.6f;
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
    real tjca=75.0f;
    jca=fcass-(fcass-jca)*exp(-dt/tjca);
    real tffp=2.5f*tff;
    ffp=fss-(fss-ffp)*exp(-dt/tffp);
    real fp=Aff*ffp+Afs*fs;
    real tfcafp=2.5f*tfcaf;
    fcafp=fcass-(fcass-fcafp)*exp(-dt/tfcafp);
    real fcap=Afcaf*fcafp+Afcas*fcas;
    real Kmn=0.002f;
    real k2n=1000.0f;
    real km2n=jca*1.0f;
    real anca=1.0f/(k2n/km2n+pow(1.0f+Kmn/cass,4.0f));
    nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
    real PhiCaL=4.0f*vffrt*(cass*exp(2.0f*vfrt)-0.341f*cao)/(exp(2.0f*vfrt)-1.0f);
    real PhiCaNa=1.0f*vffrt*(0.75f*nass*exp(1.0*vfrt)-0.75f*nao)/(exp(1.0f*vfrt)-1.0f);
    real PhiCaK=1.0f*vffrt*(0.75f*kss*exp(1.0f*vfrt)-0.75f*ko)/(exp(1.0f*vfrt)-1.0f);
    real zca=2.0f;
    real PCa=0.0001f;
    if (celltype_gpu==1)
    {
        PCa*=1.2f;
    }
    if (celltype_gpu==2)
    {
        PCa*=2.5f;
    }
    real PCap=1.1f*PCa;
    real PCaNa=0.00125f*PCa;
    real PCaK=3.574e-4*PCa;
    real PCaNap=0.00125f*PCap;
    real PCaKp=3.574e-4*PCap;
    real fICaLp=(1.0f/(1.0f+KmCaMK/CaMKa));
    real ICaL=(1.0f-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0f-nca)+jca*fcap*nca);
    real ICaNa=(1.0f-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0f-nca)+jca*fcap*nca);
    real ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0f-nca)+jca*fcap*nca);

    real xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
    real txrf=12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
    real txrs=1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
    real Axrf=1.0/(1.0+exp((v+54.81)/38.21));
    real Axrs=1.0-Axrf;
    xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
    xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
    real xr=Axrf*xrf+Axrs*xrs;
    real rkr=1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
    real GKr=0.046;
    if (celltype_gpu==1)
    {
        GKr*=1.3;
    }
    if (celltype_gpu==2)
    {
        GKr*=0.8;
    }
    real IKr = GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);

    real xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
    real txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
    xs1=xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
    real xs2ss=xs1ss;
    real txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
    xs2=xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
    real KsCa=1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
    real GKs=0.0034;
    if (celltype_gpu==1)
    {
        GKs*=1.4;
    }
    real  IKs=GKs*KsCa*xs1*xs2*(v-EKs);

    real xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
    real txk1=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
    xk1=xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
    real rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
    real GK1=0.1908;
    if (celltype_gpu==1)
    {
        GK1*=1.2;
    }
    if (celltype_gpu==2)
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
    real allo=1.0/(1.0+pow(KmCaAct/cai,2.0));
    real zna=1.0;
    real JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    real JncxCa=E2*k2-E1*k1;
    real Gncx=0.0008;
    if (celltype_gpu==1)
    {
        Gncx*=1.1;
    }
    if (celltype_gpu==2)
    {
        Gncx*=1.4;
    }

    real INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);

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
    KmCaAct=150.0e-6;
    allo=1.0/(1.0+pow(KmCaAct/cass,2.0));
    JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    JncxCa=E2*k2-E1*k1;
    real INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);

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
    real a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
    real b1=k1m*MgADP;
    real a2=k2p;
    real b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
    real a3=(k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
    real b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
    real a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
    real b4=(k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
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
    real Pnak=30;
    if (celltype_gpu==1)
    {
        Pnak*=0.9;
    }
    if (celltype_gpu==2)
    {
        Pnak*=0.7;
    }
    real INaK=Pnak*(zna*JnakNa+zk*JnakK);

    real xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
    real GKb=0.003;
    if (celltype_gpu==1)
    {
        GKb*=0.6;
    }
    real IKb=GKb*xkb*(v-EK);

    real PNab=3.75e-10;
    real INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

    real PCab=2.5e-8;
    real ICab=PCab*4.0*vffrt*(cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);

    real GpCa=0.0005;
    real IpCa=GpCa*cai/(0.0005+cai);

    //voltage()
    v+=-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+stim_current);

    CaMKb = CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
    CaMKa = CaMKb+CaMKt;
    CaMKt+=dt*(aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt);

    real JdiffNa=(nass-nai)/2.0;
    real JdiffK=(kss-ki)/2.0;
    real Jdiff=(cass-cai)/0.2;

    real bt=4.75;
    real a_rel=0.5*bt;
    real Jrel_inf=a_rel*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
    if (celltype_gpu==2)
    {
        Jrel_inf*=1.7;
    }
    real tau_rel=bt/(1.0+0.0123/cajsr);
    if (tau_rel<0.005)
    {
        tau_rel=0.005;
    }
    Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
    real btp=1.25*bt;
    real a_relp=0.5*btp;
    real Jrel_infp=a_relp*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
    if (celltype_gpu==2)
    {
        Jrel_infp*=1.7;
    }
    real tau_relp=btp/(1.0+0.0123/cajsr);
    if (tau_relp<0.005)
    {
        tau_relp=0.005;
    }
    Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp);
    real fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
    real Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;

    real Jupnp=0.004375*cai/(cai+0.00092);
    real Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);
    if (celltype_gpu==1)
    {
        Jupnp*=1.3;
        Jupp*=1.3;
    }
    real fJupp=(1.0/(1.0+KmCaMK/CaMKa));
    real Jleak=0.0039375*cansr/15.0;
    real Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

    real Jtr=(cansr-cajsr)/100.0;

    nai+=dt*(-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);
    nass+=dt*(-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);

    ki+=dt*(-(Ito+IKr+IKs+IK1+IKb+stim_current-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);
    kss+=dt*(-(ICaK)*Acap/(F*vss)-JdiffK);

    real Bcai;
    if (celltype_gpu==1)
    {
        Bcai=1.0/(1.0+1.3*cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
    }
    else
    {
        Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
    }
    cai+=dt*(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));

    real Bcass=1.0/(1.0+BSRmax*KmBSR/pow(KmBSR+cass,2.0)+BSLmax*KmBSL/pow(KmBSL+cass,2.0));
    cass+=dt*(Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff));

    cansr+=dt*(Jup-Jtr*vjsr/vnsr);

    real Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow(kmcsqn+cajsr,2.0));
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
