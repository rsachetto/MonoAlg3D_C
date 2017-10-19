#include <stddef.h>
#include <unitypes.h>
#include "../main/constants.h"
#include <stdlib.h>
#include <stdio.h>
#include "model_gpu_utils.h"

#define ENDO
#define NEQ 12


static __device__ size_t pitch;
static size_t pitch_h;

__global__ void kernel_set_model_inital_conditions(Real *sv, int num_volumes);

__global__ void solve_gpu(Real dt, Real *sv, Real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, Real *fibrosis, Real atpi);

inline __device__ void RHS_gpu(Real *sv_, Real *rDY_, Real stim_current, int threadID_, Real dt, Real fibrosis, Real atpi);


extern "C" size_t set_model_initial_conditions_gpu(Real **sv, uint32_t num_volumes) {

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(Real);

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}


extern "C" void solve_model_odes_gpu(Real dt, Real *sv, Real *stim_currents, uint32_t *cells_to_solve,
                                    uint32_t num_cells_to_solve, int num_steps, void *extra_data,
                                    size_t extra_data_bytes_size) {


    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(Real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    Real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    Real atpi = 6.8;
    Real *fibrosis_device;
    Real *fibs = NULL;

    if(extra_data) {
        atpi = ((Real*)extra_data)[0];
        fibs = ((Real*)extra_data)+1;
    }
    else {
        fibs = (Real*)calloc(num_cells_to_solve, sizeof(Real));
    }

    check_cuda_error(cudaMalloc((void **) &fibrosis_device, extra_data_bytes_size-sizeof(Real)));
    check_cuda_error(cudaMemcpy(fibrosis_device, fibs, extra_data_bytes_size-sizeof(Real), cudaMemcpyHostToDevice));

    solve_gpu<<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, fibrosis_device, atpi);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    check_cuda_error(cudaFree(fibrosis_device));

    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if(!extra_data) free(fibs);
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
                          int num_steps, Real *fibrosis,  Real atpi)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        Real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, fibrosis[threadID], atpi);

            *((Real*)((char*)sv) + sv_id) = dt*rDY[0] + *((Real*)((char*)sv) + sv_id);

            for(int i = 1; i < 12; i++) {
                *((Real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }

        }

    }
}


inline __device__ void RHS_gpu(Real *sv_, Real *rDY_, Real stim_current, int threadID_, Real dt, Real fibrosis, Real atpi) {

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
    const Real hatp = 2;             // Hill coefficient

    Real Ko   = 5.4;

    Real atpi_change = 6.8f-atpi;

    atpi = atpi + atpi_change*fibrosis;

    //Real katp = 0.306;
    const Real katp = -0.0942857142857f*atpi + 0.683142857143f; //Ref: A Comparison of Two Models of Human Ventricular Tissue: Simulated Ischaemia and Re-entry    


    const Real patp =  1.0f/(1.0f + powf((atpi/katp),hatp));
    const Real gkatp    =  0.000195f/nicholsarea;
    const Real gkbaratp =  gkatp*patp*powf((Ko/4),natp);

    const Real katp2= 1.4;
    const Real hatp2 = 2.6;
    const Real pcal = 1.0f/(1.0f + powf((katp2/atpi),hatp2));

    const Real Cao=2.0;
    const Real Nao=140.0;
    const Real Cai=0.00007;
    const Real Nai=7.67;
    const Real Ki=138.3;

//Constants
    const Real R=8314.472;
    const Real F=96485.3415f;
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
    const Real GCaL=0.2786f*pcal;
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


    const Real Ek=RTONF*(logf((Ko/Ki)));
    const Real Ena=RTONF*(logf((Nao/Nai)));
    const Real Eks=RTONF*(logf((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    const Real Eca=0.5f*RTONF*(logf((Cao/Cai)));
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
