#include "Maleckar2008.h"
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
    log_to_stdout_and_file("Using Maleckar2008 GPU model\n");

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

    if (threadID < num_volumes) {
         *((real * )((char *) sv + pitch * 0) + threadID)  = -87.169816169406;
         *((real * )((char *) sv + pitch * 1) + threadID)  = 0.001075453357;
         *((real * )((char *) sv + pitch * 2) + threadID)  = 0.990691306716;
         *((real * )((char *) sv + pitch * 3) + threadID)  = 0.993888937283;
         *((real * )((char *) sv + pitch * 4) + threadID)  = 0.000018211252;
         *((real * )((char *) sv + pitch * 5) + threadID)  = 0.979322592773;
         *((real * )((char *) sv + pitch * 6) + threadID)  = 0.001208153482;
         *((real * )((char *) sv + pitch * 7) + threadID)  = 0.000033616596;
         *((real * )((char *) sv + pitch * 8) + threadID)  = 0.004173008466;
         *((real * )((char *) sv + pitch * 9) + threadID)  = 0.015242594688;
         *((real * )((char *) sv + pitch * 10) + threadID) = 0.007074239331;
         *((real * )((char *) sv + pitch * 11) + threadID) = 0.048267587131;
         *((real * )((char *) sv + pitch * 12) + threadID) = 0.105468807033;
         *((real * )((char *) sv + pitch * 13) + threadID) = 0.00364776906;
         *((real * )((char *) sv + pitch * 14) + threadID) = 0.174403618112;
         *((real * )((char *) sv + pitch * 15) + threadID) = 0.003643592594;
         *((real * )((char *) sv + pitch * 16) + threadID) = 0.993331326442;
         *((real * )((char *) sv + pitch * 17) + threadID) = 97.505463697266;
         *((real * )((char *) sv + pitch * 18) + threadID) = 0.006679257264;
         *((real * )((char *) sv + pitch * 19) + threadID) = 11.441712311614;
         *((real * )((char *) sv + pitch * 20) + threadID) = 1.716573130685;
         *((real * )((char *) sv + pitch * 21) + threadID) = 0.226941113355;
         *((real * )((char *) sv + pitch * 22) + threadID) = 0.256752008084;
         *((real * )((char *) sv + pitch * 23) + threadID) = 104.450004990523;
         *((real * )((char *) sv + pitch * 24) + threadID) = 22.171689894953;
         *((real * )((char *) sv + pitch * 25) + threadID) = 19.864701949854;

        if(use_adpt) {
            *((real *)((char *)sv + pitch * 26) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 27) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 28) + threadID) = 0.0;    // previous dt
        }
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve,
                          uint32_t num_cells_to_solve, int num_steps) {
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

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id) {

    #define DT *((real *)((char *)sv + pitch * 30) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * 31) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * 32) + thread_id)

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
            for(int i = 0; i < NEQ; i++) {
                _k_aux__[i] = _k2__[i];
                _k2__[i] = _k1__[i];
                _k1__[i] = _k_aux__[i];
            }

            // it steps the method ahead, with euler solution
            for(int i = 0; i < NEQ; i++) {
                sv_local[i] = edos_new_euler_[i];
            }

            if(time_new + previous_dt >= final_time) {
                if((fabs(final_time - time_new) < 1.0e-5)) {
                    break;
                } else if(time_new < final_time) {
                    dt = previous_dt = final_time - time_new;
                    time_new += previous_dt;
                    break;
                } else {
                    dt = previous_dt = min_dt;
                    time_new += (final_time - time_new);
                    printf("Error: %d: %lf\n", thread_id, final_time - time_new);
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



inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int thread_id, real dt) {
    //State variables
    double var_cell__V ;  // Units: mV; Initial value: -87.169816169406
    double var_INa__xm ; // Units: dimensionless; Initial value: 0.001075453357
    double var_INa__xh ;  // Units: dimensionless; Initial value: 0.990691306716
    double var_INa__xj ; // Units: dimensionless; Initial value: 0.993888937283
    double var_ICaL__c1 ; // Units: dimensionless; Initial value: 0.000018211252
    double var_ICaL__c2 ; // Units: dimensionless; Initial value: 0.979322592773
    double var_ICaL__xi1ca ; // Units: dimensionless; Initial value: 0.001208153482
    double var_ICaL__xi1ba ; // Units: dimensionless; Initial value: 0.000033616596
    double var_ICaL__xi2ca ; // Units: dimensionless; Initial value: 0.004173008466
    double var_ICaL__xi2ba ; // Units: dimensionless; Initial value: 0.015242594688
    double var_IKr__xr ; // Units: dimensionless; Initial value: 0.007074239331
    double var_IKs__xs1 ; // Units: dimensionless; Initial value: 0.048267587131
    double var_IKs__xs2 ;  // Units: dimensionless; Initial value: 0.105468807033
    double var_Ito__xtos ;  // Units: dimensionless; Initial value: 0.00364776906
    double var_Ito__ytos ; // Units: dimensionless; Initial value: 0.174403618112
    double var_Ito__xtof ; // Units: dimensionless; Initial value: 0.003643592594
    double var_Ito__ytof ; // Units: dimensionless; Initial value: 0.993331326442
    double var_Irel__Ca_JSR ; // Units: uM; Initial value: 97.505463697266
    double var_Irel__xir ; // Units: uM_per_ms; Initial value: 0.006679257264
    double var_Na__Na_i ; // Units: mM; Initial value: 11.441712311614
    double var_Ca__Ca_dyad ; // Units: uM; Initial value: 1.716573130685
    double var_Ca__Ca_submem ; // Units: uM; Initial value: 0.226941113355
    double var_Ca__Ca_i ; // Units: uM; Initial value: 0.256752008084
    double var_Ca__Ca_NSR ; // Units: uM; Initial value: 104.450004990523
    double var_Ca__tropi ;// Units: uM; Initial value: 22.171689894953
    double var_Ca__trops ; // Units: uM; Initial value: 19.864701949854
    
    if(use_adpt) {
        var_cell__V = sv[0];
        var_INa__xm = sv[1];
        var_INa__xh = sv[2];
        var_INa__xj = sv[3];
        var_ICaL__c1 = sv[4];
        var_ICaL__c2 = sv[5];
        var_ICaL__xi1ca = sv[6];
        var_ICaL__xi1ba = sv[7];
        var_ICaL__xi2ca = sv[8];
        var_ICaL__xi2ba = sv[9];
        var_IKr__xr = sv[10];
        var_IKs__xs1 = sv[11];
        var_IKs__xs2 = sv[12];
        var_Ito__xtos = sv[13];
        var_Ito__ytos = sv[14];
        var_Ito__xtof = sv[15];
        var_Ito__ytof = sv[16];
        var_Irel__Ca_JSR = sv[17];
        var_Irel__xir = sv[18];
        var_Na__Na_i = sv[19];
        var_Ca__Ca_dyad = sv[20];
        var_Ca__Ca_submem = sv[21];
        var_Ca__Ca_i = sv[22];
        var_Ca__Ca_NSR = sv[23];
        var_Ca__tropi = sv[24];
        var_Ca__trops = sv[25];
    } else {
        var_cell__V =  *((real*)((char*)sv + pitch * 0) + thread_id);
        var_INa__xm =  *((real*)((char*)sv + pitch * 1) + thread_id);
        var_INa__xh =  *((real*)((char*)sv + pitch * 2) + thread_id);
        var_INa__xj =  *((real*)((char*)sv + pitch * 3) + thread_id);
        var_ICaL__c1 =  *((real*)((char*)sv + pitch * 4) + thread_id);
        var_ICaL__c2 =  *((real*)((char*)sv + pitch * 5) + thread_id);
        var_ICaL__xi1ca =  *((real*)((char*)sv + pitch * 6) + thread_id);
        var_ICaL__xi1ba =  *((real*)((char*)sv + pitch * 7) + thread_id);
        var_ICaL__xi2ca =  *((real*)((char*)sv + pitch * 8) + thread_id);
        var_ICaL__xi2ba =  *((real*)((char*)sv + pitch * 9) + thread_id);
        var_IKr__xr =  *((real*)((char*)sv + pitch * 10) + thread_id);
        var_IKs__xs1 =  *((real*)((char*)sv + pitch * 11) + thread_id);
        var_IKs__xs2 =  *((real*)((char*)sv + pitch * 12) + thread_id);
        var_Ito__xtos =  *((real*)((char*)sv + pitch * 13) + thread_id);
        var_Ito__ytos =  *((real*)((char*)sv + pitch * 14) + thread_id);
        var_Ito__xtof =  *((real*)((char*)sv + pitch * 15) + thread_id);
        var_Ito__ytof =  *((real*)((char*)sv + pitch * 16) + thread_id);
        var_Irel__Ca_JSR =  *((real*)((char*)sv + pitch * 17) + thread_id);
        var_Irel__xir =  *((real*)((char*)sv + pitch * 18) + thread_id);
        var_Na__Na_i =  *((real*)((char*)sv + pitch * 19) + thread_id);
        var_Ca__Ca_dyad =  *((real*)((char*)sv + pitch * 20) + thread_id);
        var_Ca__Ca_submem =  *((real*)((char*)sv + pitch * 21) + thread_id);
        var_Ca__Ca_i =  *((real*)((char*)sv + pitch * 22) + thread_id);
        var_Ca__Ca_NSR =  *((real*)((char*)sv + pitch * 23) + thread_id);
        var_Ca__tropi =  *((real*)((char*)sv + pitch * 24) + thread_id);
        var_Ca__trops =  *((real*)((char*)sv + pitch * 25) + thread_id);
    }

    #include "Maleckar2008_common.inc.c"

}

