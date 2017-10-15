#include <stddef.h>
#include <unitypes.h>
#include "../main/constants.h"
#include <stdlib.h>
#include <stdio.h>
#include "model_gpu_utils.h"
#define ENDO


static __device__ size_t pitch;
static size_t pitchh;

__global__ void kernel_set_model_inital_conditions(Real *sv, int num_volumes);

__global__ void solve_gpu(Real dt, Real *sv, Real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          Real stim_start, Real stim_dur, Real time,
                          int num_steps, int neq, void *extra_data);

__global__ void update_refinement(Real *sv, uint32_t *cells, size_t number_of_cells, int neq);

inline __device__ void RHS_gpu(Real *sv_, Real *rDY_, Real stim_current, Real time, Real stim_start, Real stim_dur, int threadID_, Real dt);


extern "C" size_t set_model_initial_conditions_gpu(Real **sv, uint32_t num_volumes, int neq) {

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(Real);

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitchh, size, (size_t )neq));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitchh, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitchh;

}


extern "C" void solve_model_ode_gpu(Real dt, Real *sv, Real *stim_currents, uint32_t *cells_to_solve,
                                    uint32_t num_cells_to_solve, Real stim_start, Real stim_dur,
                                    Real time, int num_steps, int neq, void *extra_data) {


    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(Real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    Real *stims_currents_device;
    uint32_t *cells_to_solve_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
    check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));

    solve_gpu<<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve,
            stim_start, stim_dur, time, num_steps, neq, extra_data);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    check_cuda_error(cudaFree(cells_to_solve_device));

}


__global__ void kernel_set_model_inital_conditions(Real *sv, int num_volumes)
{
    // Thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if(threadID < num_volumes) {

        *((Real *) ((char *) sv + pitch * 0) + threadID) = -86.2f;   // V;       millivolt
        *((Real *) ((char *) sv + pitch * 1) + threadID) = 0.0f; //M
        *((Real *) ((char *) sv + pitch * 2) + threadID) = 0.75; //H
        *((Real *) ((char *) sv + pitch * 3) + threadID) = 0.75; //J
        *((Real *) ((char *) sv + pitch * 4) + threadID) = 0.0f; //Xr1
        *((Real *) ((char *) sv + pitch * 5) + threadID) = 0.0f; //Xs
        *((Real *) ((char *) sv + pitch * 6) + threadID) = 1.0; //S
        *((Real *) ((char *) sv + pitch * 7) + threadID) = 1.0; //F
        *((Real *) ((char *) sv + pitch * 8) + threadID) = 1.0; //F2
        *((Real *) ((char *) sv + pitch * 9) + threadID) = 0.0; //D_INF
        *((Real *) ((char *) sv + pitch * 10) + threadID) = 0.0; //R_INF
        *((Real *) ((char *) sv + pitch * 11) + threadID) = 0.0; //Xr2_INF
    }
}


// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(Real dt, Real *sv, Real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          Real stim_start, Real stim_dur, Real time,
                          int num_steps, int neq, void *extra_data)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;
    Real t = time;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        sv_id = cells_to_solve[threadID];
        Real *rDY = (Real *)malloc(neq*sizeof(Real));

        for (int n = 0; n < num_steps; ++n) {

            RHS_gpu(sv, rDY, stim_currents[threadID], t, stim_start, stim_dur, sv_id, dt);

            *((Real*)((char*)sv) + sv_id) = dt*rDY[0] + *((Real*)((char*)sv) + sv_id);

            for(int i = 1; i < 12; i++) {
                *((Real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }

            
            t += dt;
        }
        free(rDY);

    }
}


inline __device__ void RHS_gpu(Real *sv_, Real *rDY_, Real stim_current, Real time, Real stim_start, Real stim_dur, int threadID_, Real dt) {
    const Real svolt = *((Real*)((char*)sv_ + pitch * 0) + threadID_);

    const Real sm   = *((Real*)((char*)sv_ + pitch * 1) + threadID_);
    const Real sh   = *((Real*)((char*)sv_ + pitch * 2) + threadID_);
    const Real sj   = *((Real*)((char*)sv_ + pitch * 3) + threadID_);
    const Real sxr1 = *((Real*)((char*)sv_ + pitch * 4) + threadID_);
    const Real sxs  = *((Real*)((char*)sv_ + pitch * 5) + threadID_);
    const Real ss   = *((Real*)((char*)sv_ + pitch * 6) + threadID_);
    const Real  sf  = *((Real*)((char*)sv_ + pitch * 7) + threadID_);
    const Real sf2  = *((Real*)((char*)sv_ + pitch * 8) + threadID_);
    const Real D_INF  = *((Real*)((char*)sv_ + pitch * 9) + threadID_);
    const Real Xr2_INF  = *((Real*)((char*)sv_ + pitch * 10) + threadID_);
    const Real R_INF  = *((Real*)((char*)sv_ + pitch * 11) + threadID_);

    const Real natp = 0.24;          // K dependence of ATP-sensitive K current
    const Real nicholsarea = 0.00005; // Nichol's areas (cm^2)
    //const Real atpi = 4.0;             // Intracellular ATP concentraion (mM)
    const Real hatp = 2;             // Hill coefficient
    //const Real katp = 0.306;         // Half-maximal saturation point of ATP-sensitive K current (mM)

    Real Ko   = 5.4;
    Real atpi = 6.8;
    Real katp = 0.042;

/*    if(fibrosis == 1) {
        Ko = 8.0;
        atpi = 4;
        katp = 0.306;
    }*/

//    Ko = 8.0 - 2.6*fibrosis;
//    atpi = 4.0 + 2.8*fibrosis;
 //   katp = 0.306 - 0.264*fibrosis;

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


    const Real Ek=RTONF*(log((Ko/Ki)));
    const Real Ena=RTONF*(log((Nao/Nai)));
    const Real Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    const Real Eca=0.5*RTONF*(log((Cao/Cai)));
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

    Real Ak1;
    Real Bk1;
    Real rec_iK1;
    Real rec_ipK;
    Real rec_iNaK;
    Real AM;
    Real BM;
    Real AH_1;
    Real BH_1;
    Real AH_2;
    Real BH_2;
    Real AJ_1;
    Real BJ_1;
    Real AJ_2;
    Real BJ_2;
    Real M_INF;
    Real H_INF;
    Real J_INF;
    Real TAU_M;
    Real TAU_H;
    Real TAU_J;
    Real axr1;
    Real bxr1;
    Real Xr1_INF;
    Real Xr2_INF_new;
    Real TAU_Xr1;
    Real Axs;
    Real Bxs;
    Real Xs_INF;
    Real TAU_Xs;
    Real R_INF_new;
    Real S_INF;
    Real TAU_S;
    Real Af;
    Real Bf;
    Real Cf;
    Real Af2;
    Real Bf2;
    Real Cf2;
    Real D_INF_new;
    Real TAU_F;
    Real F_INF;
    Real TAU_F2;
    Real F2_INF;
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

    Real calc_i_stim = stim_current;

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
