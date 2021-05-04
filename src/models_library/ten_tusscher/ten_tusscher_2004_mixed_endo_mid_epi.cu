#include "../../gpu_utils/gpu_utils.h"
#include <stddef.h>
#include <stdint.h>

#include "ten_tusscher_2004_mixed_endo_mid_epi.h"

__constant__ size_t pitch;
size_t pitch_h;

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    uint32_t num_volumes = solver->original_num_cells;

    #ifdef ENDO
        log_info("Using ten Tusscher 2004 ENDO GPU model\n");
    #endif

    #ifdef EPI
        log_info("Using ten Tusscher 2004 EPI GPU model\n");
    #endif

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

	// Get the mapping array
    uint32_t *mapping = NULL;
    uint32_t *mapping_device = NULL;
    if(ode_solver->ode_extra_data) 
    {
        mapping = (uint32_t*)ode_solver->ode_extra_data;
        check_cuda_error(cudaMalloc((void **)&mapping_device, ode_solver->extra_data_size));
        check_cuda_error(cudaMemcpy(mapping_device, mapping, ode_solver->extra_data_size, cudaMemcpyHostToDevice));
    }
    else 
    {
        log_error_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    solve_gpu <<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, mapping_device);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if(mapping_device) check_cuda_error(cudaFree(mapping_device));
}

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes)
{
    // Thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if(threadID < num_volumes) {

        *((real*)((char*)sv + pitch * 0) + threadID)  = INITIAL_V;   // V;       millivolt
        *((real*)((char*)sv + pitch * 1) + threadID)  = 0.f;   //M
        *((real*)((char*)sv + pitch * 2) + threadID)  = 0.75;    //H
        *((real*)((char*)sv + pitch * 3) + threadID)  = 0.75f;    //J
        *((real*)((char*)sv + pitch * 4) + threadID)  = 0.f;   //Xr1
        *((real*)((char*)sv + pitch * 5) + threadID)  = 1.f;    //Xr2
        *((real*)((char*)sv + pitch * 6) + threadID)  = 0.f;    //Xs
        *((real*)((char*)sv + pitch * 7) + threadID)  = 1.f;  //S
        *((real*)((char*)sv + pitch * 8) + threadID)  = 0.f;    //R
        *((real*)((char*)sv + pitch * 9) + threadID)  = 0.f;    //D
        *((real*)((char*)sv + pitch * 10) + threadID) = 1.f;   //F
        *((real*)((char*)sv + pitch * 11) + threadID) = 1.f; //FCa
        *((real*)((char*)sv + pitch * 12) + threadID) = 1.f;  //G
        *((real*)((char*)sv + pitch * 13) + threadID) = 0.0002;  //Cai
        *((real*)((char*)sv + pitch * 14) + threadID) = 0.2f;      //CaSR
        *((real*)((char*)sv + pitch * 15) + threadID) = 11.6f;   //Nai
        *((real*)((char*)sv + pitch * 16) + threadID) = 138.3f;    //Ki


    }
}


// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,                         
						  uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, uint32_t *mapping)
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

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, mapping[threadID]);

            for(int i = 0; i < NEQ; i++) {
                *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }
            
        }

    }
}


inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt, int mapping) {

    // State variables
    real svolt = *((real*)((char*)sv + pitch * 0) + threadID_);
    real sm    = *((real*)((char*)sv + pitch * 1) + threadID_);
    real sh    = *((real*)((char*)sv + pitch * 2) + threadID_);
    real sj    = *((real*)((char*)sv + pitch * 3) + threadID_);
    real sxr1  = *((real*)((char*)sv + pitch * 4) + threadID_);
    real sxr2  = *((real*)((char*)sv + pitch * 5) + threadID_);
    real sxs   = *((real*)((char*)sv + pitch * 6) + threadID_);
    real ss    = *((real*)((char*)sv + pitch * 7) + threadID_);
    real sr    = *((real*)((char*)sv + pitch * 8) + threadID_);
    real sd    = *((real*)((char*)sv + pitch * 9) + threadID_);
    real sf    = *((real*)((char*)sv + pitch * 10) + threadID_);
    real sfca  = *((real*)((char*)sv + pitch * 11) + threadID_);
    real sg    = *((real*)((char*)sv + pitch * 12) + threadID_);
    real Cai   = *((real*)((char*)sv + pitch * 13) + threadID_);
    real CaSR  = *((real*)((char*)sv + pitch * 14) + threadID_);
    real Nai   = *((real*)((char*)sv + pitch * 15) + threadID_);
    real Ki    = *((real*)((char*)sv + pitch * 16) + threadID_);

	#include "ten_tusscher_2004_mixed_endo_mid_epi.common.c"
    
}
