#include <stddef.h>
#include <stdint.h>

#include "stewart_aslanidi_noble_2009.h"

__constant__  size_t pitch;
__constant__  real abstol;
__constant__  real reltol;
__constant__  real max_dt;
__constant__  real min_dt;
__constant__  uint8_t use_adpt;

size_t pitch_h;

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    uint8_t use_adpt_h = (uint8_t)solver->adaptive;

    check_cuda_error(cudaMemcpyToSymbol(use_adpt, &use_adpt_h, sizeof(uint8_t)));
    log_info("Using Stewart-Aslanidi-Noble 2009 GPU model\n");

    uint32_t num_volumes = solver->original_num_cells;

    if(use_adpt_h) {
        real reltol_h = solver->rel_tol;
        real abstol_h = solver->abs_tol;
        real max_dt_h = solver->max_dt;
        real min_dt_h = solver->min_dt;

        check_cuda_error(cudaMemcpyToSymbol(reltol, &reltol_h, sizeof(real)));
        check_cuda_error(cudaMemcpyToSymbol(abstol, &abstol_h, sizeof(real)));
        check_cuda_error(cudaMemcpyToSymbol(max_dt, &max_dt_h, sizeof(real)));
        check_cuda_error(cudaMemcpyToSymbol(min_dt, &min_dt_h, sizeof(real)));
        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_info("Using Euler model to solve the ODEs\n");
    }

    // execution configuration
    const int GRID = (num_volumes + BLOCK_SIZE - 1) / BLOCK_SIZE;

    size_t size = num_volumes * sizeof(real);

    if(use_adpt_h)
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ + 3));
    else
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ));

    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

    kernel_set_model_initial_conditions<<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes);

    check_cuda_error(cudaPeekAtLastError());
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
    const int GRID = ((int)num_cells_to_solve + BLOCK_SIZE - 1) / BLOCK_SIZE;

    size_t stim_currents_size = sizeof(real) * num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t) * num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **)&stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    // the array cells to solve is passed when we are using and adapative mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **)&cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(
            cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve,
                                    num_steps);

    check_cuda_error(cudaPeekAtLastError());

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device)
        check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        *((real * )((char *) sv + pitch * 0) + threadID) = -69.1370441635924;
        *((real * )((char *) sv + pitch * 1) + threadID) = 136.781894160227;
        *((real * )((char *) sv + pitch * 2) + threadID) = 8.80420286531673;
        *((real * )((char *) sv + pitch * 3) + threadID) = 0.000101878186157052;
        *((real * )((char *) sv + pitch * 4) + threadID) = 0.0457562667986602;
        *((real * )((char *) sv + pitch * 5) + threadID) = 0.00550281999719088;
        *((real * )((char *) sv + pitch * 6) + threadID) = 0.313213286437995;
        *((real * )((char *) sv + pitch * 7) + threadID) = 0.00953708522974789;
        *((real * )((char *) sv + pitch * 8) + threadID) = 0.0417391656294997;
        *((real * )((char *) sv + pitch * 9) + threadID) = 0.190678733735145;
        *((real * )((char *) sv + pitch * 10) + threadID) = 0.238219836154029;
        *((real * )((char *) sv + pitch * 11) + threadID) = 0.000446818714055411;
        *((real * )((char *) sv + pitch * 12) + threadID) = 0.000287906256206415;
        *((real * )((char *) sv + pitch * 13) + threadID) = 0.989328560287987;
        *((real * )((char *) sv + pitch * 14) + threadID) = 0.995474890442185;
        *((real * )((char *) sv + pitch * 15) + threadID) = 0.999955429598213;
        *((real * )((char *) sv + pitch * 16) + threadID) = 0.96386101799501;
        *((real * )((char *) sv + pitch * 17) + threadID) = 0.00103618091196912;
        *((real * )((char *) sv + pitch * 18) + threadID) = 3.10836886659417;
        *((real * )((char *) sv + pitch * 19) + threadID) = 0.991580051907845;

        if(use_adpt) {
            *((real *)((char *)sv + pitch * 20) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 21) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 22) + threadID) = 0.0;    // previous dt
        }
    }
}

//Include the default solver used by all current models.
#include "../default_solvers.cu"

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real STATES[NEQ];
    if (use_adpt)
    {
        for (uint32_t i = 0; i < NEQ; i++)
            STATES[i] = sv_[i];
    }
    else
    {
        for (uint32_t i = 0; i < NEQ; i++)
            STATES[i] = *((real *)((char *)sv_ + pitch * i) + threadID_);
    }

    #include "stewart_aslanidi_noble_2009_common.inc"
}

