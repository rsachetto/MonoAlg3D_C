#include <stddef.h>
#include <unitypes.h>
#include "../main/constants.h"
#include <stdlib.h>
#include <stdio.h>
#include "model_gpu_utils.h"

static __device__ size_t pitch;
static size_t pitch_h;

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

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitch_h, size, (size_t )neq));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

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

        *((Real*)((char*)sv + pitch * 0) + threadID) = -85.23f;   // V;       millivolt
        //Set the rest of the inital conditions here!!
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

            for(int i = 1; i < 13; i++) {
                *((Real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }

            for(int i = 13; i < 19; i++) {
                *((Real *) ((char *) sv + pitch * i) + sv_id) = dt * rDY[i] + *((Real *) ((char *) sv + pitch * i) + sv_id);
            }
            
            t += dt;
        }
        free(rDY);

    }
}


inline __device__ void RHS_gpu(Real *sv_, Real *rDY_, Real stim_current, Real time, Real stim_start, Real stim_dur, int threadID_, Real dt) {


}
