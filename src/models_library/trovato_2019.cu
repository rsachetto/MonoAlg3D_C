#include "trovato_2019.h"
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
    log_to_stdout_and_file("Using Trovato_2019 GPU model\n");

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
        *((real * )((char *) sv + pitch * 0) + threadID) = -86.6814002878592;
        *((real * )((char *) sv + pitch * 1) + threadID) = 0.00505983330678751;
        *((real * )((char *) sv + pitch * 2) + threadID) = 0.000101777993438818;
        *((real * )((char *) sv + pitch * 3) + threadID) = 8.23183964616932;
        *((real * )((char *) sv + pitch * 4) + threadID) = 8.23153516580562;
        *((real * )((char *) sv + pitch * 5) + threadID) = 8.23154325237268;
        *((real * )((char *) sv + pitch * 6) + threadID) = 143.767359809132;
        *((real * )((char *) sv + pitch * 7) + threadID) = 143.767768218104;
        *((real * )((char *) sv + pitch * 8) + threadID) = 143.767769906216;
        *((real * )((char *) sv + pitch * 9) + threadID) = 4.36004404734282e-05;
        *((real * )((char *) sv + pitch * 10) + threadID) = 0.000102004317781147;
        *((real * )((char *) sv + pitch * 11) + threadID) = 1.26350902016858;
        *((real * )((char *) sv + pitch * 12) + threadID) = 1.24811940209535;
        *((real * )((char *) sv + pitch * 13) + threadID) = 1.26516959198518;
        *((real * )((char *) sv + pitch * 14) + threadID) = 0.000108240945806962;
        *((real * )((char *) sv + pitch * 15) + threadID) = 1.25045800437317e-69;
        *((real * )((char *) sv + pitch * 16) + threadID) = 0.00632661703915808;
        *((real * )((char *) sv + pitch * 17) + threadID) = 0.788611739889677;
        *((real * )((char *) sv + pitch * 18) + threadID) = 0.788545979951331;
        *((real * )((char *) sv + pitch * 19) + threadID) = 0.790474358603666;
        *((real * )((char *) sv + pitch * 20) + threadID) = 0.579693514309867;
        *((real * )((char *) sv + pitch * 21) + threadID) = 0.790947058236417;
        *((real * )((char *) sv + pitch * 22) + threadID) = 0.000241925773627233;
        *((real * )((char *) sv + pitch * 23) + threadID) = 0.463574582508218;
        *((real * )((char *) sv + pitch * 24) + threadID) = 0.240216198686475;
        *((real * )((char *) sv + pitch * 25) + threadID) = 0.000272851144435704;
        *((real * )((char *) sv + pitch * 26) + threadID) = 0.649604795721571;
        *((real * )((char *) sv + pitch * 27) + threadID) = 0.989965695822495;
        *((real * )((char *) sv + pitch * 28) + threadID) = 6.97735089296892e-09;
        *((real * )((char *) sv + pitch * 29) + threadID) = 0.999999968230738;
        *((real * )((char *) sv + pitch * 30) + threadID) = 0.926692153319136;
        *((real * )((char *) sv + pitch * 31) + threadID) = 0.99999996819573;
        *((real * )((char *) sv + pitch * 32) + threadID) = 0.999999905741936;
        *((real * )((char *) sv + pitch * 33) + threadID) = 0.999978907334662;
        *((real * )((char *) sv + pitch * 34) + threadID) = 0.999999968365903;
        *((real * )((char *) sv + pitch * 35) + threadID) = 0.999999968278239;
        *((real * )((char *) sv + pitch * 36) + threadID) = 0.00547252500964926;
        *((real * )((char *) sv + pitch * 37) + threadID) = 0.000304250912559619;
        *((real * )((char *) sv + pitch * 38) + threadID) = 0.994214357917907;
        *((real * )((char *) sv + pitch * 39) + threadID) = 0.000331691184084272;
        *((real * )((char *) sv + pitch * 40) + threadID) = 0.568716473334161;
        *((real * )((char *) sv + pitch * 41) + threadID) = 0.191165248085394;
        *((real * )((char *) sv + pitch * 42) + threadID) = 0.000222677365291219;
        *((real * )((char *) sv + pitch * 43) + threadID) = 0.233119011214908;
        *((real * )((char *) sv + pitch * 44) + threadID) = 0.997084813729909;
        *((real * )((char *) sv + pitch * 45) + threadID) = 0.466236137183558;

        if(use_adpt) 
        {
            *((real *)((char *)sv + pitch * 46) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 47) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 48) + threadID) = 0.0;    // previous dt
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

    #define DT *((real *)((char *)sv + pitch * 46) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * 47) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * 48) + thread_id)

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

    #include "trovato_common_inc.c"
}

