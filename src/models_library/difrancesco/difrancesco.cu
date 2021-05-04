#include "../../gpu_utils/gpu_utils.h"
#include <stddef.h>
#include <stdint.h>

#include "difrancesco.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    log_info("Using DiFrancesco & Noble 1985 GPU model\n");

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

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        *((real * )((char *) sv + pitch * 0) + threadID) = -87;    // V millivolt
        *((real * )((char *) sv + pitch * 1) + threadID) = 4;      // Kc millimolar
        *((real * )((char *) sv + pitch * 2) + threadID) = 140;    // Ki millimolar
        *((real * )((char *) sv + pitch * 3) + threadID) = 8;      // Nai millimolar
        *((real * )((char *) sv + pitch * 4) + threadID) = 0.2;    // y dimensionless
        *((real * )((char *) sv + pitch * 5) + threadID) = 0.01;   // x dimensionless
        *((real * )((char *) sv + pitch * 6) + threadID) = 5e-5;   // Cai millimolar
        *((real * )((char *) sv + pitch * 7) + threadID) = 1;      // s dimensionless
        *((real * )((char *) sv + pitch * 8) + threadID) = 0.01;   // m dimensionless
        *((real * )((char *) sv + pitch * 9) + threadID) = 0.8;    // h dimensionless
        *((real * )((char *) sv + pitch * 10) + threadID) = 0.005; // d dimensionless
        *((real * )((char *) sv + pitch * 11) + threadID) = 1;     // f dimensionless
        *((real * )((char *) sv + pitch * 12) + threadID) = 1;     // f2 dimensionless
        *((real * )((char *) sv + pitch * 13) + threadID) = 2;     // Ca_up millimolar
        *((real * )((char *) sv + pitch * 14) + threadID) = 1;     // Ca_rel millimolar
        *((real * )((char *) sv + pitch * 15) + threadID) = 1;     // p dimensionless

    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id);

            for(int i = 0; i < NEQ; i++) {
                *((real *) ((char *) sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *) ((char *) sv + pitch * i) + sv_id);
            }            

        }

    }
}

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_) {

    //State variables
    real STATES[16]; 
    STATES[0] = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    STATES[1] = *((real*)((char*)sv_ + pitch * 1) + threadID_);
    STATES[2] = *((real*)((char*)sv_ + pitch * 2) + threadID_);
    STATES[3] = *((real*)((char*)sv_ + pitch * 3) + threadID_);
    STATES[4] = *((real*)((char*)sv_ + pitch * 4) + threadID_);
    STATES[5] = *((real*)((char*)sv_ + pitch * 5) + threadID_);
    STATES[6] = *((real*)((char*)sv_ + pitch * 6) + threadID_);
    STATES[7] = *((real*)((char*)sv_ + pitch * 7) + threadID_);
    STATES[8] = *((real*)((char*)sv_ + pitch * 8) + threadID_);
    STATES[9] = *((real*)((char*)sv_ + pitch * 9) + threadID_);
    STATES[10] = *((real*)((char*)sv_ + pitch * 10) + threadID_);
    STATES[11] = *((real*)((char*)sv_ + pitch * 11) + threadID_);
    STATES[12] = *((real*)((char*)sv_ + pitch * 12) + threadID_);
    STATES[13] = *((real*)((char*)sv_ + pitch * 13) + threadID_);
    STATES[14] = *((real*)((char*)sv_ + pitch * 14) + threadID_);
    STATES[15] = *((real*)((char*)sv_ + pitch * 15) + threadID_);

    // Constants
    real CONSTANTS[50]; 
    CONSTANTS[0] = 8314.472;
    CONSTANTS[1] = 310;
    CONSTANTS[2] = 96485.3415;
    CONSTANTS[3] = 0.075;
    CONSTANTS[4] = stim_current;
    CONSTANTS[5] = 3;
    CONSTANTS[6] = 3;
    CONSTANTS[7] = 45;
    CONSTANTS[8] = 140;
    CONSTANTS[9] = 1e-5;
    CONSTANTS[10] = 180;
    CONSTANTS[11] = 920;
    CONSTANTS[12] = 210;
    CONSTANTS[13] = 10;
    CONSTANTS[14] = 0.0005;
    CONSTANTS[15] = 0.28;
    CONSTANTS[16] = 0.18;
    CONSTANTS[17] = 0.02;
    CONSTANTS[18] = 2;
    CONSTANTS[19] = 125;
    CONSTANTS[20] = 1;
    CONSTANTS[21] = 40;
    CONSTANTS[22] = 3;
    CONSTANTS[23] = 0.02;
    CONSTANTS[24] = 0.001;
    CONSTANTS[25] = 0.5;
    CONSTANTS[26] = 750;
    CONSTANTS[27] = 1e-5;
    CONSTANTS[28] = 15;
    CONSTANTS[29] = 0.0001;
    CONSTANTS[30] = 0.0001;
    CONSTANTS[31] = 5;
    CONSTANTS[32] = 0.001;
    CONSTANTS[33] = 0.05;
    CONSTANTS[34] = 2;
    CONSTANTS[35] = 0.1;
    CONSTANTS[36] = 5;
    CONSTANTS[37] = 0.001;
    CONSTANTS[38] = 0.025;
    CONSTANTS[39] = 2;
    CONSTANTS[40] = 0.05;
    CONSTANTS[41] = 2;
    CONSTANTS[42] = 0.00157;
    CONSTANTS[43] = 4;
    CONSTANTS[44] = 0.7;
    CONSTANTS[45] = ( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2];
    CONSTANTS[46] =  3.14159*pow(CONSTANTS[33], 2.00000)*CONSTANTS[34];
    CONSTANTS[47] =  CONSTANTS[46]*(1.00000 - CONSTANTS[35]);
    CONSTANTS[48] =  CONSTANTS[47]*0.0500000;
    CONSTANTS[49] = CONSTANTS[47]*0.0200000;

    // Algebraics
    real ALGEBRAIC[46];
    ALGEBRAIC[8] = ( STATES[6]*CONSTANTS[31])/CONSTANTS[32];
    ALGEBRAIC[2] = ( 0.500000*exp( 0.0826000*(STATES[0]+50.0000)))/(1.00000+exp( 0.0570000*(STATES[0]+50.0000)));
    ALGEBRAIC[12] = ( 1.30000*exp( - 0.0600000*(STATES[0]+20.0000)))/(1.00000+exp( - 0.0400000*(STATES[0]+20.0000)));
    ALGEBRAIC[3] =  0.0330000*exp(- STATES[0]/17.0000);
    ALGEBRAIC[13] = 33.0000/(1.00000+exp(- (STATES[0]+10.0000)/8.00000));
    ALGEBRAIC[5] =  20.0000*exp( - 0.125000*(STATES[0]+75.0000));
    ALGEBRAIC[15] = 2000.00/( 320.000*exp( - 0.100000*(STATES[0]+75.0000))+1.00000);
    ALGEBRAIC[9] = ( 0.625000*(STATES[0]+34.0000))/(exp((STATES[0]+34.0000)/4.00000) - 1.00000);
    ALGEBRAIC[18] = 5.00000/(1.00000+exp(( - 1.00000*(STATES[0]+34.0000))/4.00000));
    ALGEBRAIC[1] =  0.0500000*exp( - 0.0670000*((STATES[0]+52.0000) - 10.0000));
    ALGEBRAIC[11] = (STATES[0]+52.0000) - 10.0000;
    ALGEBRAIC[19] = (fabs(ALGEBRAIC[11])<CONSTANTS[9] ? 2.50000 : ( 1.00000*ALGEBRAIC[11])/(1.00000 - exp( - 0.200000*ALGEBRAIC[11])));
    ALGEBRAIC[4] = STATES[0]+41.0000;
    ALGEBRAIC[14] = (fabs(ALGEBRAIC[4])<CONSTANTS[27] ? 2000.00 : ( 200.000*ALGEBRAIC[4])/(1.00000 - exp( - 0.100000*ALGEBRAIC[4])));
    ALGEBRAIC[20] =  8000.00*exp( - 0.0560000*(STATES[0]+66.0000));
    ALGEBRAIC[6] = (STATES[0]+24.0000) - 5.00000;
    ALGEBRAIC[16] = (fabs(ALGEBRAIC[6])<CONSTANTS[29] ? 120.000 : ( 30.0000*ALGEBRAIC[6])/(1.00000 - exp(( - 1.00000*ALGEBRAIC[6])/4.00000)));
    ALGEBRAIC[21] = (fabs(ALGEBRAIC[6])<CONSTANTS[29] ? 120.000 : ( 12.0000*ALGEBRAIC[6])/(exp(ALGEBRAIC[6]/10.0000) - 1.00000));
    ALGEBRAIC[7] = STATES[0]+34.0000;
    ALGEBRAIC[17] = (fabs(ALGEBRAIC[7])<CONSTANTS[30] ? 25.0000 : ( 6.25000*ALGEBRAIC[7])/(exp(ALGEBRAIC[7]/4.00000) - 1.00000));
    ALGEBRAIC[22] = 50.0000/(1.00000+exp(( - 1.00000*(STATES[0]+34.0000))/4.00000));
    ALGEBRAIC[0] =  CONSTANTS[45]*log(CONSTANTS[8]/STATES[3]);
    ALGEBRAIC[30] =  CONSTANTS[16]*(STATES[0] - ALGEBRAIC[0]);
    ALGEBRAIC[33] = ( (( CONSTANTS[19]*STATES[1])/(CONSTANTS[20]+STATES[1]))*STATES[3])/(CONSTANTS[21]+STATES[3]);
    ALGEBRAIC[34] = ( CONSTANTS[23]*( exp(( CONSTANTS[25]*(CONSTANTS[22] - 2.00000)*STATES[0])/CONSTANTS[45])*pow(STATES[3], CONSTANTS[22])*CONSTANTS[18] -  exp(( (CONSTANTS[25] - 1.00000)*(CONSTANTS[22] - 2.00000)*STATES[0])/CONSTANTS[45])*pow(CONSTANTS[8], CONSTANTS[22])*STATES[6]))/( (1.00000+ CONSTANTS[24]*( STATES[6]*pow(CONSTANTS[8], CONSTANTS[22])+ CONSTANTS[18]*pow(STATES[3], CONSTANTS[22])))*(1.00000+STATES[6]/0.00690000));
    ALGEBRAIC[35] =  CONSTANTS[45]*log((CONSTANTS[8]+ 0.120000*STATES[1])/(STATES[3]+ 0.120000*STATES[2]));
    ALGEBRAIC[36] =  CONSTANTS[26]*pow(STATES[8], 3.00000)*STATES[9]*(STATES[0] - ALGEBRAIC[35]);
    ALGEBRAIC[23] =  (( STATES[4]*STATES[1])/(STATES[1]+CONSTANTS[7]))*CONSTANTS[5]*(STATES[0] - ALGEBRAIC[0]);
    ALGEBRAIC[40] =  (( 0.0100000*CONSTANTS[28]*(STATES[0] - 50.0000))/( CONSTANTS[45]*(1.00000 - exp(( - 1.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))))*( STATES[3]*exp(50.0000/CONSTANTS[45]) -  CONSTANTS[8]*exp(( - 1.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))*STATES[10]*STATES[11]*STATES[12];
    ALGEBRAIC[39] =  (( 2.00000*1.00000*CONSTANTS[47]*CONSTANTS[2])/( 1.00000*CONSTANTS[38]*CONSTANTS[36]))*STATES[6]*(CONSTANTS[36] - STATES[13]);
    ALGEBRAIC[41] =  (( 2.00000*1.00000*CONSTANTS[49]*CONSTANTS[2])/( 1.00000*CONSTANTS[39]))*STATES[15]*(STATES[13] - STATES[14]);
    ALGEBRAIC[26] = ( CONSTANTS[10]*(STATES[2] -  STATES[1]*exp(- STATES[0]/CONSTANTS[45])))/140.000;
    ALGEBRAIC[27] =  STATES[5]*ALGEBRAIC[26];
    ALGEBRAIC[10] =  CONSTANTS[45]*log(STATES[1]/STATES[2]);
    ALGEBRAIC[28] = ( (( CONSTANTS[11]*STATES[1])/(STATES[1]+CONSTANTS[12]))*(STATES[0] - ALGEBRAIC[10]))/(1.00000+exp(( ((STATES[0]+10.0000) - ALGEBRAIC[10])*2.00000)/CONSTANTS[45]));
    ALGEBRAIC[29] =  (( (( STATES[7]*CONSTANTS[15]*(0.200000+STATES[1]/(CONSTANTS[13]+STATES[1]))*STATES[6])/(CONSTANTS[14]+STATES[6]))*(STATES[0]+10.0000))/(1.00000 - exp( - 0.200000*(STATES[0]+10.0000))))*( STATES[2]*exp(( 0.500000*STATES[0])/CONSTANTS[45]) -  STATES[1]*exp(( - 0.500000*STATES[0])/CONSTANTS[45]));
    ALGEBRAIC[24] =  (( STATES[4]*STATES[1])/(STATES[1]+CONSTANTS[7]))*CONSTANTS[6]*(STATES[0] - ALGEBRAIC[10]);
    ALGEBRAIC[38] =  (( 0.0100000*CONSTANTS[28]*(STATES[0] - 50.0000))/( CONSTANTS[45]*(1.00000 - exp(( - 1.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))))*( STATES[2]*exp(50.0000/CONSTANTS[45]) -  STATES[1]*exp(( - 1.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))*STATES[10]*STATES[11]*STATES[12];
    ALGEBRAIC[42] = (ALGEBRAIC[28]+ALGEBRAIC[27]+ALGEBRAIC[24]+ALGEBRAIC[38]+ALGEBRAIC[29]) -  2.00000*ALGEBRAIC[33];
    ALGEBRAIC[25] = ALGEBRAIC[23]+ALGEBRAIC[24];
    ALGEBRAIC[31] =  0.500000*CONSTANTS[45]*log(CONSTANTS[18]/STATES[6]);
    ALGEBRAIC[32] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[31]);
    ALGEBRAIC[37] =  (( 4.00000*CONSTANTS[28]*(STATES[0] - 50.0000))/( CONSTANTS[45]*(1.00000 - exp(( - 1.00000*(STATES[0] - 50.0000)*2.00000)/CONSTANTS[45]))))*( STATES[6]*exp(100.000/CONSTANTS[45]) -  CONSTANTS[18]*exp(( - 2.00000*(STATES[0] - 50.0000))/CONSTANTS[45]))*STATES[10]*STATES[11]*STATES[12];
    ALGEBRAIC[43] = ALGEBRAIC[37]+ALGEBRAIC[38]+ALGEBRAIC[40];
    ALGEBRAIC[44] = ( (( 2.00000*1.00000*CONSTANTS[49]*CONSTANTS[2])/( 1.00000*CONSTANTS[40]))*STATES[14]*pow(STATES[6], CONSTANTS[41]))/(pow(STATES[6], CONSTANTS[41])+pow(CONSTANTS[37], CONSTANTS[41]));

    real RATES[16];
    RATES[0] = (- (ALGEBRAIC[25]+ALGEBRAIC[27]+ALGEBRAIC[28]+ALGEBRAIC[29]+ALGEBRAIC[30]+ALGEBRAIC[32]+ALGEBRAIC[33]+ALGEBRAIC[34]+ALGEBRAIC[36]+ALGEBRAIC[43]+CONSTANTS[4])/CONSTANTS[3]) * 1.0E-03;
    RATES[1] =  (- CONSTANTS[44]*(STATES[1] - CONSTANTS[43])+( 1.00000*ALGEBRAIC[42])/( 1.00000*CONSTANTS[42]*CONSTANTS[2])) * 1.0E-03;
    RATES[2] = ( ( - 1.00000*ALGEBRAIC[42])/( 1.00000*CONSTANTS[47]*CONSTANTS[2]) ) * 1.0E-03;
    RATES[3] = ( ( - 1.00000*(ALGEBRAIC[36]+ALGEBRAIC[30]+ALGEBRAIC[23]+ALGEBRAIC[40]+ ALGEBRAIC[33]*3.00000+( ALGEBRAIC[34]*CONSTANTS[22])/(CONSTANTS[22] - 2.00000)))/( 1.00000*CONSTANTS[47]*CONSTANTS[2])) * 1.0E-03;
    RATES[4] =  (ALGEBRAIC[1]*(1.00000 - STATES[4]) -  ALGEBRAIC[19]*STATES[4]) * 1.0E-03;
    RATES[5] =  (ALGEBRAIC[2]*(1.00000 - STATES[5]) -  ALGEBRAIC[12]*STATES[5]) * 1.0E-03;
    RATES[6] = (( - 1.00000*((((ALGEBRAIC[37]+ALGEBRAIC[32]) - ( 2.00000*ALGEBRAIC[34])/(CONSTANTS[22] - 2.00000)) - ALGEBRAIC[44])+ALGEBRAIC[39]))/( 2.00000*1.00000*CONSTANTS[47]*CONSTANTS[2])) * 1.0E-03;
    RATES[7] =  (ALGEBRAIC[3]*(1.00000 - STATES[7]) -  ALGEBRAIC[13]*STATES[7]) * 1.0E-03;
    RATES[8] =  (ALGEBRAIC[14]*(1.00000 - STATES[8]) -  ALGEBRAIC[20]*STATES[8]) * 1.0E-03;
    RATES[9] =  (ALGEBRAIC[5]*(1.00000 - STATES[9]) -  ALGEBRAIC[15]*STATES[9]) * 1.0E-03;
    RATES[10] =  (ALGEBRAIC[16]*(1.00000 - STATES[10]) -  ALGEBRAIC[21]*STATES[10]) * 1.0E-03;
    RATES[11] =  (ALGEBRAIC[17]*(1.00000 - STATES[11]) -  ALGEBRAIC[22]*STATES[11]) * 1.0E-03;
    RATES[12] = (CONSTANTS[31] -  STATES[12]*(CONSTANTS[31]+ALGEBRAIC[8])) * 1.0E-03;
    RATES[13] = (( 1.00000*(ALGEBRAIC[39] - ALGEBRAIC[41]))/( 2.00000*1.00000*CONSTANTS[48]*CONSTANTS[2])) * 1.0E-03;
    RATES[14] = (( 1.00000*(ALGEBRAIC[41] - ALGEBRAIC[44]))/( 2.00000*1.00000*CONSTANTS[49]*CONSTANTS[2])) * 1.0E-03;
    RATES[15] = (ALGEBRAIC[9]*(1.00000 - STATES[15]) - ALGEBRAIC[18]*STATES[15]) * 1.0E-03;

    // Rates
    rDY_[0] = RATES[0];
    rDY_[1] = RATES[1];
    rDY_[2] = RATES[2];
    rDY_[3] = RATES[3];
    rDY_[4] = RATES[4];
    rDY_[5] = RATES[5];
    rDY_[6] = RATES[6];
    rDY_[7] = RATES[7];
    rDY_[8] = RATES[8];
    rDY_[9] = RATES[9];
    rDY_[10] = RATES[10];
    rDY_[11] = RATES[11];
    rDY_[12] = RATES[12];
    rDY_[13] = RATES[13];
    rDY_[14] = RATES[14];
    rDY_[15] = RATES[15];

    
}

