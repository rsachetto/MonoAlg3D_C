#include <stddef.h>
#include "../monodomain/constants.h"
#include "model_gpu_utils.h"

#include "ten_tusscher_3_RS.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    print_to_stdout_and_file("Using ten Tusscher 3 GPU model\n");

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}

extern "C" SOLVE_MODEL_ODES_GPU(solve_model_odes_gpu) {

    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    // Default values for a healthy cell ///////////
    real atpi = 6.8f;
    real Ko = 5.4f;
    real Ki_mult = 1.0f;
    real acidosis = 0.0;
    real K1_mult = 1.0f;
    ////////////////////////////////////

    real *fibrosis_device;
    real *fibs = NULL;

    if(extra_data) {
        atpi = ((real*)extra_data)[0]; //value
        Ko = ((real*)extra_data)[1]; //value
        Ki_mult = ((real*)extra_data)[2]; //value
        K1_mult = ((real*)extra_data)[3]; //value
        acidosis = ((real*)extra_data)[4]; //value
        fibs = ((real*)extra_data) + 5; //pointer

        extra_data_bytes_size = extra_data_bytes_size-5*sizeof(real);
    }
    else {
        extra_data_bytes_size = num_cells_to_solve*sizeof(real);
        fibs = (real*)calloc(num_cells_to_solve, sizeof(real));
    }

    check_cuda_error(cudaMalloc((void **) &fibrosis_device, extra_data_bytes_size));
    check_cuda_error(cudaMemcpy(fibrosis_device, fibs, extra_data_bytes_size, cudaMemcpyHostToDevice));

    solve_gpu<<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, fibrosis_device, atpi, Ko, Ki_mult, K1_mult, acidosis);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    check_cuda_error(cudaFree(fibrosis_device));

    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if(!extra_data) free(fibs);
}


__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes)
{
    // Thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if(threadID < num_volumes) {

        *((real *) ((char *) sv + pitch * 0) + threadID) = -86.2f;   // V;       millivolt
        *((real *) ((char *) sv + pitch * 1) + threadID) = 0.0f; //M
        *((real *) ((char *) sv + pitch * 2) + threadID) = 0.75; //H
        *((real *) ((char *) sv + pitch * 3) + threadID) = 0.75; //J
        *((real *) ((char *) sv + pitch * 4) + threadID) = 0.0f; //Xr1
        *((real *) ((char *) sv + pitch * 5) + threadID) = 0.0f; //Xs
        *((real *) ((char *) sv + pitch * 6) + threadID) = 1.0; //S
        *((real *) ((char *) sv + pitch * 7) + threadID) = 1.0; //F
        *((real *) ((char *) sv + pitch * 8) + threadID) = 1.0; //F2
        *((real *) ((char *) sv + pitch * 9) + threadID) = 0.0; //D_INF
        *((real *) ((char *) sv + pitch * 10) + threadID) = 0.0; //R_INF
        *((real *) ((char *) sv + pitch * 11) + threadID) = 0.0; //Xr2_INF
    }
}


// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, real *fibrosis,  real atpi,
                          real Ko, real Ki_multiplicator, real K1_multiplicator, real acidosis)
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

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, fibrosis[threadID], atpi, Ko, Ki_multiplicator, K1_multiplicator, acidosis);

            *((real*)((char*)sv) + sv_id) = dt*rDY[0] + *((real*)((char*)sv) + sv_id);

            for(int i = 1; i < 12; i++) {
                *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }

        }

    }
}


inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt, real fibrosis, real atpi, real Ko, real Ki_multiplicator, real K1_multiplicator, real acidosis) {
    //fibrosis = 0 means that the cell is fibrotic, 1 is not fibrotic. Anything between 0 and 1 means border zone
    const real svolt = *((real*)((char*)sv_ + pitch * 0) + threadID_);

    //printf("%lf, %lf, %lf, %lf, %lf, %lf\n", atpi, Ko, Ki_multiplicator, acidosis, fibrosis, K1_multiplicator);

    real svolt_acid = svolt;

    if( (fibrosis == 0.0f) && (acidosis == 1.0f) ) {
        //These values are from In Electrophysiologic effects of acute myocardial ischemia: a theoretical
        //study of altered cell excitability and action potential duration
        svolt_acid = svolt - 3.4f;
    }

    const real sm   = *((real*)((char*)sv_ + pitch * 1) + threadID_);
    const real sh   = *((real*)((char*)sv_ + pitch * 2) + threadID_);
    const real sj   = *((real*)((char*)sv_ + pitch * 3) + threadID_);
    const real sxr1 = *((real*)((char*)sv_ + pitch * 4) + threadID_);
    const real sxs  = *((real*)((char*)sv_ + pitch * 5) + threadID_);
    const real ss   = *((real*)((char*)sv_ + pitch * 6) + threadID_);
    const real sf  = *((real*)((char*)sv_ + pitch * 7) + threadID_);
    const real sf2  = *((real*)((char*)sv_ + pitch * 8) + threadID_);
    const real D_INF  = *((real*)((char*)sv_ + pitch * 9) + threadID_);
    const real R_INF  = *((real*)((char*)sv_ + pitch * 10) + threadID_);    
    const real Xr2_INF  = *((real*)((char*)sv_ + pitch * 11) + threadID_);


    #include "ten_tusscher_3_common.inc"

}
