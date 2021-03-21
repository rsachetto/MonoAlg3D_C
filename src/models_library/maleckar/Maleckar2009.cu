#include "Maleckar2009.h"
#include <stddef.h>
#include <stdint.h>

__constant__ size_t pitch;
__constant__ real abstol;
__constant__ real reltol;
__constant__ real max_dt;
__constant__ real min_dt;
__constant__ uint8_t use_adpt;
size_t pitch_h;

#define sv(i) *((real *)((char *)sv + pitch * (i)) + thread_id)

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    uint8_t use_adpt_h = (uint8_t)solver->adaptive;

    check_cuda_error(cudaMemcpyToSymbol(use_adpt, &use_adpt_h, sizeof(uint8_t)));
    log_info("Using Maleckar2009 GPU model\n");

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
    uint32_t *cells_to_solve = ode_solver->cells_to_solve;
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
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error(cudaPeekAtLastError());

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device)
        check_cuda_error(cudaFree(cells_to_solve_device));
}

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes) {

	int thread_id = blockDim.x * blockIdx.x + threadIdx.x;

    if(thread_id < num_volumes) {

        sv(0) = -73.941851;
        sv(1) = 0.003325;
        sv(2) = 0.875262;
        sv(3) = 0.870692;
        sv(4) = 0.000014;
        sv(5) = 0.998578;
        sv(6) = 0.998561;
        sv(7) = 0.001098;
        sv(8) = 0.948202;
        sv(9) = 0.000371;
        sv(10) = 0.966869;
        sv(11) = 0.004661;
        sv(12) = 0.000054;
        sv(13) = 8.488527;
        sv(14) = 6.5e-5;
        sv(15) = 129.502075;
        sv(16) = 7.1e-5;
        sv(17) = 0.026604;
        sv(18) = 0.012843;
        sv(19) = 0.190077;
        sv(20) = 0.714719;
        sv(21) = 1.38222;
        sv(22) = 130.019282;
        sv(23) = 1.814418;
        sv(24) = 5.588239;
        sv(25) = 0.630471;
        sv(26) = 0.646226;
        sv(27) = 0.43071;
        sv(28) = 0.45453;
        sv(29) = 0.002665;

        if(use_adpt) {
            sv(NEQ) = min_dt; // dt
            sv(NEQ+1) = 0.0;    // time_new
            sv(NEQ+2) = 0.0;    // previous dt
        }
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, uint32_t num_cells_to_solve, int num_steps) {
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
                    *((real *)((char *)sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
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
    // State variables
    real var_membrane__V;                                          // Units: millivolt; Initial value: -73.941851
    real var_sodium_current_m_gate__m;                             // Units: dimensionless; Initial value: 0.003325
    real var_sodium_current_h1_gate__h1;                           // Units: dimensionless; Initial value: 0.875262
    real var_sodium_current_h2_gate__h2;                           // Units: dimensionless; Initial value: 0.870692
    real var_L_type_Ca_channel_d_L_gate__d_L;                      // Units: dimensionless; Initial value: 0.000014
    real var_L_type_Ca_channel_f_L1_gate__f_L1;                    // Units: dimensionless; Initial value: 0.998578
    real var_L_type_Ca_channel_f_L2_gate__f_L2;                    // Units: dimensionless; Initial value: 0.998561
    real var_Ca_independent_transient_outward_K_current_r_gate__r; // Units: dimensionless; Initial value: 0.001098
    real var_Ca_independent_transient_outward_K_current_s_gate__s; // Units: dimensionless; Initial value: 0.948202
    real var_ultra_rapid_K_current_aur_gate__a_ur;                 // Units: dimensionless; Initial value: 0.000371
    real var_ultra_rapid_K_current_iur_gate__i_ur;                 // Units: dimensionless; Initial value: 0.966869
    real var_delayed_rectifier_K_currents_n_gate__n;               // Units: dimensionless; Initial value: 0.004661
    real var_delayed_rectifier_K_currents_pa_gate__pa;             // Units: dimensionless; Initial value: 0.000054
    real var_intracellular_ion_concentrations__Na_i;               // Units: millimolar; Initial value: 8.488527
    real var_intracellular_ion_concentrations__Ca_i;               // Units: millimolar; Initial value: 6.5e-5
    real var_intracellular_ion_concentrations__K_i;                // Units: millimolar; Initial value: 129.502075
    real var_intracellular_ion_concentrations__Ca_d;               // Units: millimolar; Initial value: 7.1e-5
    real var_intracellular_Ca_buffering__O_C;                      // Units: dimensionless; Initial value: 0.026604
    real var_intracellular_Ca_buffering__O_TC;                     // Units: dimensionless; Initial value: 0.012843
    real var_intracellular_Ca_buffering__O_TMgC;                   // Units: dimensionless; Initial value: 0.190077
    real var_intracellular_Ca_buffering__O_TMgMg;                  // Units: dimensionless; Initial value: 0.714719
    real var_cleft_space_ion_concentrations__Na_c;                 // Units: millimolar; Initial value: 130.019282
    real var_cleft_space_ion_concentrations__Ca_c;                 // Units: millimolar; Initial value: 1.814418
    real var_cleft_space_ion_concentrations__K_c;                  // Units: millimolar; Initial value: 5.588239
    real var_Ca_handling_by_the_SR__Ca_rel;                        // Units: millimolar; Initial value: 0.630471
    real var_Ca_handling_by_the_SR__Ca_up;                         // Units: millimolar; Initial value: 0.646226
    real var_Ca_handling_by_the_SR__O_Calse;                       // Units: dimensionless; Initial value: 0.43071
    real var_Ca_handling_by_the_SR__F1;                            // Units: dimensionless; Initial value: 0.45453
    real var_Ca_handling_by_the_SR__F2;                            // Units: dimensionless; Initial value: 0.002665

    if(use_adpt) {
        var_membrane__V = sv[0];                                          // Units: millivolt; Initial value: -73.941851
        var_sodium_current_m_gate__m = sv[1];                             // Units: dimensionless; Initial value: 0.003325
        var_sodium_current_h1_gate__h1 = sv[2];                           // Units: dimensionless; Initial value: 0.875262
        var_sodium_current_h2_gate__h2 = sv[3];                           // Units: dimensionless; Initial value: 0.870692
        var_L_type_Ca_channel_d_L_gate__d_L = sv[4];                      // Units: dimensionless; Initial value: 0.000014
        var_L_type_Ca_channel_f_L1_gate__f_L1 = sv[5];                    // Units: dimensionless; Initial value: 0.998578
        var_L_type_Ca_channel_f_L2_gate__f_L2 = sv[6];                    // Units: dimensionless; Initial value: 0.998561
        var_Ca_independent_transient_outward_K_current_r_gate__r = sv[7]; // Units: dimensionless; Initial value: 0.001098
        var_Ca_independent_transient_outward_K_current_s_gate__s = sv[8]; // Units: dimensionless; Initial value: 0.948202
        var_ultra_rapid_K_current_aur_gate__a_ur = sv[9];                 // Units: dimensionless; Initial value: 0.000371
        var_ultra_rapid_K_current_iur_gate__i_ur = sv[10];                // Units: dimensionless; Initial value: 0.966869
        var_delayed_rectifier_K_currents_n_gate__n = sv[11];              // Units: dimensionless; Initial value: 0.004661
        var_delayed_rectifier_K_currents_pa_gate__pa = sv[12];            // Units: dimensionless; Initial value: 0.000054
        var_intracellular_ion_concentrations__Na_i = sv[13];              // Units: millimolar; Initial value: 8.488527
        var_intracellular_ion_concentrations__Ca_i = sv[14];              // Units: millimolar; Initial value: 6.5e-5
        var_intracellular_ion_concentrations__K_i = sv[15];               // Units: millimolar; Initial value: 129.502075
        var_intracellular_ion_concentrations__Ca_d = sv[16];              // Units: millimolar; Initial value: 7.1e-5
        var_intracellular_Ca_buffering__O_C = sv[17];                     // Units: dimensionless; Initial value: 0.026604
        var_intracellular_Ca_buffering__O_TC = sv[18];                    // Units: dimensionless; Initial value: 0.012843
        var_intracellular_Ca_buffering__O_TMgC = sv[19];                  // Units: dimensionless; Initial value: 0.190077
        var_intracellular_Ca_buffering__O_TMgMg = sv[20];                 // Units: dimensionless; Initial value: 0.714719
        var_cleft_space_ion_concentrations__Na_c = sv[22];                // Units: millimolar; Initial value: 130.019282
        var_cleft_space_ion_concentrations__Ca_c = sv[23];                // Units: millimolar; Initial value: 1.814418
        var_cleft_space_ion_concentrations__K_c = sv[24];                 // Units: millimolar; Initial value: 5.588239
        var_Ca_handling_by_the_SR__Ca_rel = sv[25];                       // Units: millimolar; Initial value: 0.630471
        var_Ca_handling_by_the_SR__Ca_up = sv[26];                        // Units: millimolar; Initial value: 0.646226
        var_Ca_handling_by_the_SR__O_Calse = sv[27];                      // Units: dimensionless; Initial value: 0.43071
        var_Ca_handling_by_the_SR__F1 = sv[28];                           // Units: dimensionless; Initial value: 0.45453
        var_Ca_handling_by_the_SR__F2 = sv[29];                           // Units: dimensionless; Initial value: 0.002665

    } else {
        var_membrane__V = sv(0);                                          // Units: millivolt; Initial value: -73.941851
        var_sodium_current_m_gate__m = sv(1);                             // Units: dimensionless; Initial value: 0.003325
        var_sodium_current_h1_gate__h1 = sv(2);                           // Units: dimensionless; Initial value: 0.875262
        var_sodium_current_h2_gate__h2 = sv(3);                           // Units: dimensionless; Initial value: 0.870692
        var_L_type_Ca_channel_d_L_gate__d_L = sv(4);                      // Units: dimensionless; Initial value: 0.000014
        var_L_type_Ca_channel_f_L1_gate__f_L1 = sv(5);                    // Units: dimensionless; Initial value: 0.998578
        var_L_type_Ca_channel_f_L2_gate__f_L2 = sv(6);                    // Units: dimensionless; Initial value: 0.998561
        var_Ca_independent_transient_outward_K_current_r_gate__r = sv(7); // Units: dimensionless; Initial value: 0.001098
        var_Ca_independent_transient_outward_K_current_s_gate__s = sv(8); // Units: dimensionless; Initial value: 0.948202
        var_ultra_rapid_K_current_aur_gate__a_ur = sv(9);                 // Units: dimensionless; Initial value: 0.000371
        var_ultra_rapid_K_current_iur_gate__i_ur = sv(10);                // Units: dimensionless; Initial value: 0.966869
        var_delayed_rectifier_K_currents_n_gate__n = sv(11);              // Units: dimensionless; Initial value: 0.004661
        var_delayed_rectifier_K_currents_pa_gate__pa = sv(12);            // Units: dimensionless; Initial value: 0.000054
        var_intracellular_ion_concentrations__Na_i = sv(13);              // Units: millimolar; Initial value: 8.488527
        var_intracellular_ion_concentrations__Ca_i = sv(14);              // Units: millimolar; Initial value: 6.5e-5
        var_intracellular_ion_concentrations__K_i = sv(15);               // Units: millimolar; Initial value: 129.502075
        var_intracellular_ion_concentrations__Ca_d = sv(16);              // Units: millimolar; Initial value: 7.1e-5
        var_intracellular_Ca_buffering__O_C = sv(17);                     // Units: dimensionless; Initial value: 0.026604
        var_intracellular_Ca_buffering__O_TC = sv(18);                    // Units: dimensionless; Initial value: 0.012843
        var_intracellular_Ca_buffering__O_TMgC = sv(19);                  // Units: dimensionless; Initial value: 0.190077
        var_intracellular_Ca_buffering__O_TMgMg = sv(20);                 // Units: dimensionless; Initial value: 0.714719
        var_cleft_space_ion_concentrations__Na_c = sv(22);                // Units: millimolar; Initial value: 130.019282
        var_cleft_space_ion_concentrations__Ca_c = sv(23);                // Units: millimolar; Initial value: 1.814418
        var_cleft_space_ion_concentrations__K_c = sv(24);                 // Units: millimolar; Initial value: 5.588239
        var_Ca_handling_by_the_SR__Ca_rel = sv(25);                       // Units: millimolar; Initial value: 0.630471
        var_Ca_handling_by_the_SR__Ca_up = sv(26);                        // Units: millimolar; Initial value: 0.646226
        var_Ca_handling_by_the_SR__O_Calse = sv(27);                      // Units: dimensionless; Initial value: 0.43071
        var_Ca_handling_by_the_SR__F1 = sv(28);                           // Units: dimensionless; Initial value: 0.45453
        var_Ca_handling_by_the_SR__F2 = sv(29);                           // Units: dimensionless; Initial value: 0.002665
    }

	#include "Maleckar2009_common.inc.c"
}
