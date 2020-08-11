#include "new_torord_fkatp_endo.h"
#include <stddef.h>
#include <stdint.h>

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
    log_to_stdout_and_file("Using New_ToRORd_fkatp_endo GPU model\n");

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
        log_to_stdout_and_file("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_to_stdout_and_file("Using Euler model to solve the ODEs\n");
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

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) 
    {
        real STATES[43];
        STATES[0] = -88.7638;
        STATES[1] = 0.0111;
        STATES[2] = 7.0305e-5;
        STATES[3] = 12.1025;
        STATES[4] = 12.1029;
        STATES[5] = 142.3002;
        STATES[6] = 142.3002;
        STATES[7] = 1.5211;
        STATES[8] = 1.5214;
        STATES[9] = 8.1583e-05;
        STATES[10] = 8.0572e-4;
        STATES[11] = 0.8286;
        STATES[12] = 0.8284;
        STATES[13] = 0.6707;
        STATES[14] = 0.8281;
        STATES[15] = 1.629e-4;
        STATES[16] = 0.5255;
        STATES[17] = 0.2872;
        STATES[18] = 9.5098e-4;
        STATES[19] = 0.9996;
        STATES[20] = 0.5936;
        STATES[21] = 4.8454e-4;
        STATES[22] = 0.9996;
        STATES[23] = 0.6538;
        STATES[24] = 8.1084e-9;
        STATES[25] = 1.0;
        STATES[26] = 0.939;
        STATES[27] = 1.0;
        STATES[28] = 0.9999;
        STATES[29] = 1.0;
        STATES[30] = 1.0;
        STATES[31] = 1.0;
        STATES[32] = 6.6462e-4;
        STATES[33] = 0.0012;
        STATES[34] = 7.0344e-4;
        STATES[35] = 8.5109e-4;
        STATES[36] = 0.9981;
        STATES[37] = 1.3289e-5;
        STATES[38] = 3.7585e-4;
        STATES[39] = 0.248;
        STATES[40] = 1.7707e-4;
        STATES[41] = 1.6129e-22;
        STATES[42] = 1.2475e-20;

        for (int i = 0; i < NEQ; i++)
            *((real * )((char *) sv + pitch * i) + threadID) = STATES[i];

        if(use_adpt) 
        {
            *((real *)((char *)sv + pitch * 43) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 44) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 45) + threadID) = 0.0;    // previous dt
        }
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                            uint32_t *cells_to_solve, uint32_t num_cells_to_solve, int num_steps) 
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        if(!use_adpt) {
            real rDY[NEQ];

            for(int n = 0; n < num_steps; ++n) {

                RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt);

                for(int i = 0; i < NEQ; i++) {
                    *((real *)((char *)sv + pitch * i) + sv_id) =
                        dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
                }
            }
        } else {
            solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], cur_time + max_dt, sv_id);
        }
    }
}

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id) 
{

    #define DT *((real *)((char *)sv + pitch * 43) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * 44) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * 45) + thread_id)

    real rDY[NEQ];

    real _tolerances_[NEQ];
    real _aux_tol = 0.0;
    real dt = DT;
    real time_new = TIME_NEW;
    real previous_dt = PREVIOUS_DT;

    real edos_old_aux_[NEQ];
    real edos_new_euler_[NEQ];
    real _k1__[NEQ];
    real _k2__[NEQ];
    real _k_aux__[NEQ];
    real sv_local[NEQ];

    const real _beta_safety_ = 0.8;

    const real __tiny_ = powf(abstol, 2.0f);

    // dt = ((time_new + dt) > final_time) ? (final_time - time_new) : dt;
    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    //#pragma unroll
    for(int i = 0; i < NEQ; i++) {
        sv_local[i] = *((real *)((char *)sv + pitch * i) + thread_id);
    }

    RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt);
    time_new += dt;

    //#pragma unroll
    for(int i = 0; i < NEQ; i++) {
        _k1__[i] = rDY[i];
    }

    int count = 0;

    int count_limit = (final_time - time_new) / min_dt;

    int aux_count_limit = count_limit + 2000000;

    if(aux_count_limit > 0) {
        count_limit = aux_count_limit;
    }

    while(1) {

        for(int i = 0; i < NEQ; i++) {
            // stores the old variables in a vector
            edos_old_aux_[i] = sv_local[i];
            // //computes euler method
            edos_new_euler_[i] = _k1__[i] * dt + edos_old_aux_[i];
            // steps ahead to compute the rk2 method
            sv_local[i] = edos_new_euler_[i];
        }

        time_new += dt;

        RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt);
        time_new -= dt; // step back

        real greatestError = 0.0, auxError = 0.0;
        //#pragma unroll
        for(int i = 0; i < NEQ; i++) {

            // stores the new evaluation
            _k2__[i] = rDY[i];
            _aux_tol = fabs(edos_new_euler_[i]) * reltol;
            _tolerances_[i] = (abstol > _aux_tol) ? abstol : _aux_tol;

            // finds the greatest error between  the steps
            auxError = fabs(((dt / 2.0) * (_k1__[i] - _k2__[i])) / _tolerances_[i]);

            greatestError = (auxError > greatestError) ? auxError : greatestError;
        }

        /// adapt the time step
        greatestError += __tiny_;
        previous_dt = dt;
        /// adapt the time step
        dt = _beta_safety_ * dt * sqrt(1.0f / greatestError);

        if(time_new + dt > final_time) {
            dt = final_time - time_new;
        }

        // it doesn't accept the solution
        if(count < count_limit && (greatestError >= 1.0f)) {
            // restore the old values to do it again
            for(int i = 0; i < NEQ; i++) {
                sv_local[i] = edos_old_aux_[i];
            }
            count++;
            // throw the results away and compute again
        } else {
            count = 0;

            // if(greatestError >=1.0) {
            //    printf("Thread //d,accepting solution with error > //lf \n", threadID, greatestError);
            //}

            // it accepts the solutions
            // int aux = (dt > max_step && max_step != 0);
            // dt = (aux) ? max_step : dt;

            if(dt < min_dt) {
                dt = min_dt;
            }

            else if(dt > max_dt && max_dt != 0) {
                dt = max_dt;
            }

            if(time_new + dt > final_time) {
                dt = final_time - time_new;
            }

            // change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração
            //#pragma unroll
            for(int i = 0; i < NEQ; i++) {
                _k_aux__[i] = _k2__[i];
                _k2__[i] = _k1__[i];
                _k1__[i] = _k_aux__[i];
            }

            // it steps the method ahead, with euler solution
            //#pragma unroll
            for(int i = 0; i < NEQ; i++) {
                sv_local[i] = edos_new_euler_[i];
            }

            // verifica se o incremento para a próxima iteração ultrapassa o tempo de salvar, q neste caso é o tempo
            // final
            if(time_new + previous_dt >= final_time) {
                // se são iguais, ja foi calculada a iteração no ultimo passo de tempo e deve-se para o laço
                // nao usar igualdade - usar esta conta, pode-se mudar a tolerância
                // printf("//d: //lf\n", threadID, fabs(final_time - time_new));
                if((fabs(final_time - time_new) < 1.0e-5)) {
                    break;
                } else if(time_new < final_time) {
                    dt = previous_dt = final_time - time_new;
                    time_new += previous_dt;
                    break;
                } else {
                    dt = previous_dt = min_dt;
                    time_new += (final_time - time_new);
                    printf("Nao era pra chegar aqui: %d: %lf\n", thread_id, final_time - time_new);
                    break;
                }
            } else {
                time_new += previous_dt;
            }
        }
    }

    //#pragma unroll
    for(int i = 0; i < NEQ; i++) {
        *((real *)((char *)sv + pitch * i) + thread_id) = sv_local[i];
    }

    DT = dt;
    TIME_NEW = time_new;
    PREVIOUS_DT = previous_dt;
}

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt) 
{
    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real STATES[NEQ];
    if (use_adpt)
    {
        for (uint32_t i = 0; i < NEQ; i++)
            STATES[i] = sv[i];
    }
    else
    {
        for (uint32_t i = 0; i < NEQ; i++)
            STATES[i] = *((real *)((char *)sv + pitch * i) + threadID_);
    }

    #include "new_torord_fkatp_endo_common.inc.c"
}

