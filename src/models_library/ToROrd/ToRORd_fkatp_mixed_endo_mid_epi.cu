#include "ToRORd_fkatp_mixed_endo_mid_epi.h"
#include <stddef.h>
#include <stdint.h>

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    // Default initial conditions (endocardium cell)
    if (threadID < num_volumes) {
        for (int i = 0; i < NEQ; i++) {
            *((real * )((char *) sv + pitch * 0) + threadID) = -88.7638;
            *((real * )((char *) sv + pitch * 1) + threadID) = 0.0111;
            *((real * )((char *) sv + pitch * 2) + threadID)= 7.0305e-5;
            *((real * )((char *) sv + pitch * 3) + threadID)= 12.1025;
            *((real * )((char *) sv + pitch * 4) + threadID)= 12.1029;
            *((real * )((char *) sv + pitch * 5) + threadID)= 142.3002;
            *((real * )((char *) sv + pitch * 6) + threadID)= 142.3002;
            *((real * )((char *) sv + pitch * 7) + threadID)= 1.5211;
            *((real * )((char *) sv + pitch * 8) + threadID)= 1.5214;
            *((real * )((char *) sv + pitch * 9) + threadID)= 8.1583e-05;
            *((real * )((char *) sv + pitch * 10) + threadID) = 8.0572e-4;
            *((real * )((char *) sv + pitch * 11) + threadID) = 0.8286;
            *((real * )((char *) sv + pitch * 12) + threadID) = 0.8284;
            *((real * )((char *) sv + pitch * 13) + threadID) = 0.6707;
            *((real * )((char *) sv + pitch * 14) + threadID) = 0.8281;
            *((real * )((char *) sv + pitch * 15) + threadID) = 1.629e-4;
            *((real * )((char *) sv + pitch * 16) + threadID) = 0.5255;
            *((real * )((char *) sv + pitch * 17) + threadID) = 0.2872;
            *((real * )((char *) sv + pitch * 18) + threadID) = 9.5098e-4;
            *((real * )((char *) sv + pitch * 19) + threadID) = 0.9996;
            *((real * )((char *) sv + pitch * 20) + threadID) = 0.5936;
            *((real * )((char *) sv + pitch * 21) + threadID) = 4.8454e-4;
            *((real * )((char *) sv + pitch * 22) + threadID) = 0.9996;
            *((real * )((char *) sv + pitch * 23) + threadID) = 0.6538;
            *((real * )((char *) sv + pitch * 24) + threadID) = 8.1084e-9;
            *((real * )((char *) sv + pitch * 25) + threadID) = 1.0;
            *((real * )((char *) sv + pitch * 26) + threadID) = 0.939;
            *((real * )((char *) sv + pitch * 27) + threadID) = 1.0;
            *((real * )((char *) sv + pitch * 28) + threadID) = 0.9999;
            *((real * )((char *) sv + pitch * 29) + threadID) = 1.0;
            *((real * )((char *) sv + pitch * 30) + threadID) = 1.0;
            *((real * )((char *) sv + pitch * 31) + threadID) = 1.0;
            *((real * )((char *) sv + pitch * 32) + threadID) = 6.6462e-4;
            *((real * )((char *) sv + pitch * 33) + threadID) = 0.0012;
            *((real * )((char *) sv + pitch * 34) + threadID) = 7.0344e-4;
            *((real * )((char *) sv + pitch * 35) + threadID) = 8.5109e-4;
            *((real * )((char *) sv + pitch * 36) + threadID) = 0.9981;
            *((real * )((char *) sv + pitch * 37) + threadID) = 1.3289e-5;
            *((real * )((char *) sv + pitch * 38) + threadID) = 3.7585e-4;
            *((real * )((char *) sv + pitch * 39) + threadID) = 0.248;
            *((real * )((char *) sv + pitch * 40) + threadID) = 1.7707e-4;
            *((real * )((char *) sv + pitch * 41) + threadID) = 1.6129e-22;
            *((real * )((char *) sv + pitch * 42) + threadID) = 1.2475e-20;
        }
            
        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * 43) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 44) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 45) + threadID) = 0.0;    // previous dt
        }
    }
}

__global__ void kernel_set_model_initial_conditions_endo_mid_epi(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt,\
                                                real *initial_endo, real *initial_epi, real *initial_mid, real *mapping) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {
        for (int i = 0; i < NEQ; i++) {
            if (mapping[threadID] == ENDO)
                *((real * )((char *) sv + pitch * i) + threadID) = initial_endo[i];
            else if (mapping[threadID] == EPI)
                *((real * )((char *) sv + pitch * i) + threadID) = initial_epi[i];
            else
                *((real * )((char *) sv + pitch * i) + threadID) = initial_mid[i];
        }
            
        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * 43) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 44) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 45) + threadID) = 0.0;    // previous dt
        }
    }
}

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    size_t pitch_h;

    uint8_t use_adpt_dt = (uint8_t)solver->adaptive;

    log_info("Using GPU model implemented in %s\n", __FILE__);

    uint32_t num_volumes = solver->original_num_cells;

    if(use_adpt_dt) {
        log_info("Using Adaptive timestep to solve the ODEs\n");
    } else {
        log_info("Using Fixed timestep to solve the ODEs\n");
    }

    // execution configuration
    const int GRID = (num_volumes + BLOCK_SIZE - 1) / BLOCK_SIZE;

    size_t size = num_volumes * sizeof(real);

    if(use_adpt_dt)
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ + 3));
    else
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ));

    // Get initial condition from extra_data
    real *initial_conditions_endo = NULL;
    real *initial_conditions_epi = NULL;
    real *initial_conditions_mid = NULL;
    real *mapping = NULL;
    real *initial_conditions_endo_device = NULL;
    real *initial_conditions_epi_device = NULL;
    real *initial_conditions_mid_device = NULL;
    real *mapping_device = NULL;

    if(solver->ode_extra_data) {
        initial_conditions_endo = (real *)solver->ode_extra_data;
        initial_conditions_epi = (real *)solver->ode_extra_data+NEQ;
        initial_conditions_mid = (real *)solver->ode_extra_data+NEQ+NEQ;
        mapping = (real *)solver->ode_extra_data+NEQ+NEQ+NEQ;
        check_cuda_error(cudaMalloc((void **)&initial_conditions_endo_device, sizeof(real)*NEQ));
        check_cuda_error(cudaMemcpy(initial_conditions_endo_device, initial_conditions_endo, sizeof(real)*NEQ, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&initial_conditions_epi_device, sizeof(real)*NEQ));
        check_cuda_error(cudaMemcpy(initial_conditions_epi_device, initial_conditions_epi, sizeof(real)*NEQ, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&initial_conditions_mid_device, sizeof(real)*NEQ));
        check_cuda_error(cudaMemcpy(initial_conditions_mid_device, initial_conditions_mid, sizeof(real)*NEQ, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&mapping_device, sizeof(real)*num_volumes));
        check_cuda_error(cudaMemcpy(mapping_device, mapping, sizeof(real)*num_volumes, cudaMemcpyHostToDevice));
    }
    else {
        log_info("[INFO] You should supply a mask function to tag the cells when using this mixed model!\n");
        log_info("[INFO] Considering all cells ENDO!\n");
    }

    if (solver->ode_extra_data) {
        kernel_set_model_initial_conditions_endo_mid_epi<<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes, pitch_h, use_adpt_dt, solver->min_dt,\
                                                            initial_conditions_endo_device, initial_conditions_epi_device, initial_conditions_mid_device,\
                                                            mapping_device);
    }
    else {
        kernel_set_model_initial_conditions<<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes, pitch_h, use_adpt_dt, solver->min_dt);
    }
    

    check_cuda_error(cudaPeekAtLastError());
    cudaDeviceSynchronize();

    check_cuda_error(cudaFree(initial_conditions_endo_device));
    check_cuda_error(cudaFree(initial_conditions_epi_device));
    check_cuda_error(cudaFree(initial_conditions_mid_device));
    check_cuda_error(cudaFree(mapping_device));

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

    // the array cells to solve is passed when we are using and adaptive mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **)&cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    // Get the mapping array if exists
    uint32_t num_volumes = ode_solver->original_num_cells;
    real *mapping = NULL;
    real *mapping_device = NULL;
    if(ode_solver->ode_extra_data) {
        mapping = (real *)ode_solver->ode_extra_data+NEQ+NEQ+NEQ;
        check_cuda_error(cudaMalloc((void **)&mapping_device, sizeof(real)*num_volumes));
        check_cuda_error(cudaMemcpy(mapping_device, mapping, sizeof(real)*num_volumes, cudaMemcpyHostToDevice));
    }

    // Transmurality mapping defined on 'extra_data' function
    if (ode_solver->ode_extra_data) {
        solve_endo_mid_epi_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, mapping_device,\
                                    num_cells_to_solve, num_steps, ode_solver->pitch, ode_solver->adaptive, ode_solver->abs_tol, ode_solver->rel_tol, ode_solver->max_dt);
    }
    // No transmurality: all cells ENDO
    else {
        solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device,\
                                    num_cells_to_solve, num_steps, ode_solver->pitch, ode_solver->adaptive, ode_solver->abs_tol, ode_solver->rel_tol, ode_solver->max_dt);
    }
    
    check_cuda_error(cudaPeekAtLastError());

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if (mapping_device) check_cuda_error(cudaFree(mapping_device));
}

__global__ void solve_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve,\
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt) {
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
            real a[NEQ], b[NEQ];

            // TODO: Remove this boolean array and write the full algebraics instead ...
            // -------------------------------------------------------------------------------------------
            // MODEL SPECIFIC:
            // set the variables which are non-linear and hodkin-huxley type
            const real TOLERANCE = 1e-8;
            bool is_rush_larsen[NEQ];
            for (int i = 0; i < NEQ; i++) {
                is_rush_larsen[i] = ((i >= 10 && i <= 31) || (i >= 39 && i <= 42)) ? true : false;        
            }
            // -------------------------------------------------------------------------------------------

            for(int n = 0; n < num_steps; ++n) {

                RHS_RL_gpu(a, b, sv, rDY, stim_currents[threadID], 0.0, sv_id, dt, pitch, false);

                // Solve variables based on its type:
                //  Non-linear = Euler
                //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
                for (int i = 0; i < NEQ; i++) {
                    real rY = *((real *)((char *)sv + pitch * i) + sv_id);
                    if (is_rush_larsen[i]) {
                        if (abs(a[i]) < TOLERANCE) { 
                            *((real *)((char *)sv + pitch * i) + sv_id) = rY + dt * ( rY * a[i] + b[i]);
                        } 
                        else {
                            real aux = b[i] / a[i];
                            *((real *)((char *)sv + pitch * i) + sv_id) = exp(a[i] * dt)*(rY + aux) - aux;
                        }
                    }
                    else {
                        *((real *)((char *)sv + pitch * i) + sv_id) = rY + dt * rDY[i];
                    }
                }
            }
        } else {
            //solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], 0.0, cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
            solve_rush_larsen_gpu_adpt(sv, stim_currents[threadID], 0.0, cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
        }
    }
}

__global__ void solve_endo_mid_epi_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, real *mapping,\
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt) {
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
            real a[NEQ], b[NEQ];

            // TODO: Remove this boolean array and write the full algebraics instead ...
            // -------------------------------------------------------------------------------------------
            // MODEL SPECIFIC:
            // set the variables which are non-linear and hodkin-huxley type
            const real TOLERANCE = 1e-8;
            bool is_rush_larsen[NEQ];
            for (int i = 0; i < NEQ; i++) {
                is_rush_larsen[i] = ((i >= 10 && i <= 31) || (i >= 39 && i <= 42)) ? true : false;        
            }
            // -------------------------------------------------------------------------------------------

            for(int n = 0; n < num_steps; ++n) {

                RHS_RL_gpu(a, b, sv, rDY, stim_currents[threadID], mapping[threadID], sv_id, dt, pitch, false);

                // Solve variables based on its type:
                //  Non-linear = Euler
                //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
                for (int i = 0; i < NEQ; i++) {
                    real rY = *((real *)((char *)sv + pitch * i) + sv_id);
                    if (is_rush_larsen[i]) {
                        if (abs(a[i]) < TOLERANCE) { 
                            *((real *)((char *)sv + pitch * i) + sv_id) = rY + dt * ( rY * a[i] + b[i]);
                        } 
                        else {
                            real aux = b[i] / a[i];
                            *((real *)((char *)sv + pitch * i) + sv_id) = exp(a[i] * dt)*(rY + aux) - aux;
                        }
                    }
                    else {
                        *((real *)((char *)sv + pitch * i) + sv_id) = rY + dt * rDY[i];
                    }
                }
            }
        } else {
            //solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], mapping[threadID], cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
            solve_rush_larsen_gpu_adpt(sv, stim_currents[threadID], mapping[threadID], cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
        }
    }
}

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt) {

    #define DT *((real *)((char *)sv + pitch * (NEQ)) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * (NEQ+1)) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * (NEQ+2)) + thread_id)

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

    const real __tiny_ = pow(abstol, 2.0);

    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    for(int i = 0; i < NEQ; i++) {
        sv_local[i] = *((real *)((char *)sv + pitch * i) + thread_id);
    }

    RHS_gpu(sv_local, rDY, stim_curr, mapping, thread_id, dt, pitch, true);
    time_new += dt;

    for(int i = 0; i < NEQ; i++) {
        _k1__[i] = rDY[i];
    }

	while(1) {

		for(int i = 0; i < NEQ; i++) {
			// stores the old variables in a vector
			edos_old_aux_[i] = sv_local[i];
			// computes euler method
			edos_new_euler_[i] = _k1__[i] * dt + edos_old_aux_[i];
			// steps ahead to compute the rk2 method
			sv_local[i] = edos_new_euler_[i];
		}

		time_new += dt;

		RHS_gpu(sv_local, rDY, stim_curr, mapping, thread_id, dt, pitch, true);
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

		if(dt < min_dt) {
			dt = min_dt;
		}
		else if(dt > max_dt) {
			dt = max_dt;
		}

		if(time_new + dt > final_time) {
			dt = final_time - time_new;
		}

		// it doesn't accept the solution or accept and risk a NaN
		if(greatestError >= 1.0f && dt > min_dt) {
			// restore the old values to do it again
			for(int i = 0; i < NEQ; i++) {
				sv_local[i] = edos_old_aux_[i];
			}
		
		} else {
			for(int i = 0; i < NEQ; i++) {
				_k_aux__[i] = _k2__[i];
				_k2__[i] = _k1__[i];
				_k1__[i] = _k_aux__[i];
			}

			for(int i = 0; i < NEQ; i++) {
				sv_local[i] = edos_new_euler_[i];
			}

			if(time_new + previous_dt >= final_time) {
				if(final_time == time_new) {
					break;
				} else if(time_new < final_time) {
					dt = previous_dt = final_time - time_new;
					time_new += previous_dt;
					break;
				} 	
			} else {
				time_new += previous_dt;
			}
		}
	}

    for(int i = 0; i < NEQ; i++) {
        *((real *)((char *)sv + pitch * i) + thread_id) = sv_local[i];
    }

    DT = dt;
    TIME_NEW = time_new;
    PREVIOUS_DT = previous_dt;
}

inline __device__ void solve_rush_larsen_gpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt) {

    #define DT *((real *)((char *)sv + pitch * (NEQ)) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * (NEQ+1)) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * (NEQ+2)) + thread_id)

    // TODO: Remove this boolean array and write the full algebraics instead ...
    // -------------------------------------------------------------------------------------------
    // MODEL SPECIFIC:
    // set the variables which are non-linear and hodkin-huxley type
    const real TOLERANCE = 1e-8;
    bool is_rush_larsen[NEQ];
    for (int i = 0; i < NEQ; i++) {
        is_rush_larsen[i] = ((i >= 10 && i <= 31) || (i >= 39 && i <= 42)) ? true : false;        
    }
    // -------------------------------------------------------------------------------------------

    real rDY[NEQ], a_[NEQ], b_[NEQ], a_new[NEQ], b_new[NEQ];

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

    const real __tiny_ = pow(abstol, 2.0);

    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    for(int i = 0; i < NEQ; i++) {
        sv_local[i] = *((real *)((char *)sv + pitch * i) + thread_id);
    }

    RHS_RL_gpu(a_, b_, sv_local, rDY, stim_curr, mapping, thread_id, dt, pitch, true);
    time_new += dt;

    for(int i = 0; i < NEQ; i++) {
        _k1__[i] = rDY[i];
    }

	while(1) {

		for(int i = 0; i < NEQ; i++) {
			// stores the old variables in a vector
			edos_old_aux_[i] = sv_local[i];
			// computes euler method
			if (is_rush_larsen[i])
                edos_new_euler_[i] = (a_[i] < TOLERANCE) ? edos_old_aux_[i] + (edos_old_aux_[i] * a_[i] + b_[i])*(dt) : \
                                                  exp(a_[i]*(dt))*(edos_old_aux_[i] + (b_[i] / a_[i])) - (b_[i] / a_[i]);
            else
                edos_new_euler_[i] = _k1__[i] * dt + edos_old_aux_[i];
			// steps ahead to compute the rk2 method
			sv_local[i] = edos_new_euler_[i];
		}

		time_new += dt;

		RHS_RL_gpu(a_new, b_new, sv_local, rDY, stim_curr, mapping, thread_id, dt, pitch, true);
		time_new -= dt; // step back

		real greatestError = 0.0, auxError = 0.0;
		
		for(int i = 0; i < NEQ; i++) {
            // stores the new evaluation
            _k2__[i] = rDY[i];

            // Check if the variable is Hodkin-Huxley or Non-linear type
            if (is_rush_larsen[i]) {
                real as = (a_[i] + a_new[i]) * 0.5;
                real bs = (b_[i] + b_new[i]) * 0.5;
                real y_2nd_order = (fabs(as) < TOLERANCE) ? edos_old_aux_[i] + (dt) * (edos_old_aux_[i]*as + bs) : \
                                                       exp(as*(dt))*(edos_old_aux_[i] + (bs/as)) - (bs/as);
                auxError = (fabs(y_2nd_order) < TOLERANCE) ? fabs(edos_new_euler_[i] - TOLERANCE) : \
                                                        fabs( (y_2nd_order - edos_new_euler_[i])/(y_2nd_order) );
                // finds the greatest error between  the steps
                greatestError = (auxError > greatestError) ? auxError : greatestError;
            }
            else {
                real f = (_k1__[i] + _k2__[i]) * 0.5;
                real y_2nd_order = edos_old_aux_[i] + (dt) * f;
                auxError = (fabs(y_2nd_order) < TOLERANCE) ? fabs(edos_new_euler_[i] - TOLERANCE) : \
                                                        fabs( (y_2nd_order - edos_new_euler_[i])/(y_2nd_order) );
                // finds the greatest error between  the steps
                greatestError = (auxError > greatestError) ? auxError : greatestError;
            }
		}

		/// adapt the time step
		greatestError += __tiny_;
		previous_dt = dt;

		/// adapt the time step
		//dt = _beta_safety_ * dt * sqrt(1.0f / greatestError);        // Sachetto`s formula
        dt = dt * sqrt(0.5 * reltol / greatestError);                  // Jhonny`s formula

		if(dt < min_dt) {
			dt = min_dt;
		}
		else if(dt > max_dt) {
			dt = max_dt;
		}

		if(time_new + dt > final_time) {
			dt = final_time - time_new;
		}

		// it doesn't accept the solution or accept and risk a NaN
		if(greatestError >= 1.0f && dt > min_dt) {
			// restore the old values to do it again
			for(int i = 0; i < NEQ; i++) {
				sv_local[i] = edos_old_aux_[i];
			}
		
		} else {
			for(int i = 0; i < NEQ; i++) {
				_k_aux__[i] = _k2__[i];
				_k2__[i] = _k1__[i];
				_k1__[i] = _k_aux__[i];

                _k_aux__[i] = a_[i];
                a_[i] = a_new[i];
                a_new[i] = _k_aux__[i];

                _k_aux__[i] = b_[i];
                b_[i] = b_new[i];
                b_new[i] = _k_aux__[i];
			}

			for(int i = 0; i < NEQ; i++) {
				sv_local[i] = edos_new_euler_[i];
			}

			if(time_new + previous_dt >= final_time) {
				if(final_time == time_new) {
					break;
				} else if(time_new < final_time) {
					dt = previous_dt = final_time - time_new;
					time_new += previous_dt;
					break;
				} 	
			} else {
				time_new += previous_dt;
			}
		}
	}

    for(int i = 0; i < NEQ; i++) {
        *((real *)((char *)sv + pitch * i) + thread_id) = sv_local[i];
    }

    DT = dt;
    TIME_NEW = time_new;
    PREVIOUS_DT = previous_dt;
}

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, real mapping, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {
    
    // Get the celltype for the current cell
    real celltype = mapping;
    
    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real v;
    real CaMKt;
    real cass;
    real nai;
    real nass;
    real ki;
    real kss;
    real cansr;
    real cajsr;
    real cai;
    real m;
    real h;
    real j;
    real hp;
    real jp;
    real mL;
    real hL;
    real hLp;
    real a;
    real iF;
    real iS;
    real ap;
    real iFp;
    real iSp;
    real d;
    real ff;
    real fs;
    real fcaf;
    real fcas;
    real jca;
    real ffp;
    real fcafp;
    real nca_ss;
    real nca_i;
    real C1;
    real C2;
    real C3;
    real I;
    real O;
    real xs1;
    real xs2;
    real Jrel_np;
    real Jrel_p;

    if (use_adpt_dt) {
        v = sv[0];
        CaMKt = sv[1];
        cass = sv[2];
        nai = sv[3];
        nass = sv[4];
        ki = sv[5];
        kss = sv[6];
        cansr = sv[7];
        cajsr = sv[8];
        cai = sv[9];
        m = sv[10];
        h = sv[11];
        j = sv[12];
        hp = sv[13];
        jp = sv[14];
        mL = sv[15];
        hL = sv[16];
        hLp = sv[17];
        a = sv[18];
        iF = sv[19];
        iS = sv[20];
        ap = sv[21];
        iFp = sv[22];
        iSp = sv[23];
        d = sv[24];
        ff = sv[25];
        fs = sv[26];
        fcaf = sv[27];
        fcas = sv[28];
        jca = sv[29];
        ffp = sv[30];
        fcafp = sv[31];
        nca_ss = sv[32];
        nca_i = sv[33];
        C1 = sv[34];
        C2 = sv[35];
        C3 = sv[36];
        I = sv[37];
        O = sv[38];
        xs1 = sv[39];
        xs2 = sv[40];
        Jrel_np = sv[41];
        Jrel_p = sv[42];
    } else {
        v = *((real *)((char *)sv + pitch * 0) + threadID_);
        CaMKt = *((real *)((char *)sv + pitch * 1) + threadID_);
        cass = *((real *)((char *)sv + pitch * 2) + threadID_);
        nai = *((real *)((char *)sv + pitch * 3) + threadID_);
        nass = *((real *)((char *)sv + pitch * 4) + threadID_);
        ki = *((real *)((char *)sv + pitch * 5) + threadID_);
        kss = *((real *)((char *)sv + pitch * 6) + threadID_);
        cansr = *((real *)((char *)sv + pitch * 7) + threadID_);
        cajsr = *((real *)((char *)sv + pitch * 8) + threadID_);
        cai = *((real *)((char *)sv + pitch * 9) + threadID_);
        m = *((real *)((char *)sv + pitch * 10) + threadID_);
        h = *((real *)((char *)sv + pitch * 11) + threadID_);
        j = *((real *)((char *)sv + pitch * 12) + threadID_);
        hp = *((real *)((char *)sv + pitch * 13) + threadID_);
        jp = *((real *)((char *)sv + pitch * 14) + threadID_);
        mL = *((real *)((char *)sv + pitch * 15) + threadID_);
        hL = *((real *)((char *)sv + pitch * 16) + threadID_);
        hLp = *((real *)((char *)sv + pitch * 17) + threadID_);
        a = *((real *)((char *)sv + pitch * 18) + threadID_);
        iF = *((real *)((char *)sv + pitch * 19) + threadID_);
        iS = *((real *)((char *)sv + pitch * 20) + threadID_);
        ap = *((real *)((char *)sv + pitch * 21) + threadID_);
        iFp = *((real *)((char *)sv + pitch * 22) + threadID_);
        iSp = *((real *)((char *)sv + pitch * 23) + threadID_);
        d = *((real *)((char *)sv + pitch * 24) + threadID_);
        ff = *((real *)((char *)sv + pitch * 25) + threadID_);
        fs = *((real *)((char *)sv + pitch * 26) + threadID_);
        fcaf = *((real *)((char *)sv + pitch * 27) + threadID_);
        fcas = *((real *)((char *)sv + pitch * 28) + threadID_);
        jca = *((real *)((char *)sv + pitch * 29) + threadID_);
        ffp = *((real *)((char *)sv + pitch * 30) + threadID_);
        fcafp = *((real *)((char *)sv + pitch * 31) + threadID_);
        nca_ss = *((real *)((char *)sv + pitch * 32) + threadID_);
        nca_i = *((real *)((char *)sv + pitch * 33) + threadID_);
        C1 = *((real *)((char *)sv + pitch * 34) + threadID_);
        C2 = *((real *)((char *)sv + pitch * 35) + threadID_);
        C3 = *((real *)((char *)sv + pitch * 36) + threadID_);
        I = *((real *)((char *)sv + pitch * 37) + threadID_);
        O = *((real *)((char *)sv + pitch * 38) + threadID_);
        xs1 = *((real *)((char *)sv + pitch * 39) + threadID_);
        xs2 = *((real *)((char *)sv + pitch * 40) + threadID_);
        Jrel_np = *((real *)((char *)sv + pitch * 41) + threadID_);
        Jrel_p = *((real *)((char *)sv + pitch * 42) + threadID_);
    }

    #include "ToRORd_fkatp_mixed_endo_mid_epi.common.c"
}

inline __device__ void RHS_RL_gpu(real *a_, real *b_, real *sv, real *rDY_, real stim_current, real mapping, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {
    
    // Get the celltype for the current cell
    real celltype = mapping;
    
    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real v;
    real CaMKt;
    real cass;
    real nai;
    real nass;
    real ki;
    real kss;
    real cansr;
    real cajsr;
    real cai;
    real m;
    real h;
    real j;
    real hp;
    real jp;
    real mL;
    real hL;
    real hLp;
    real a;
    real iF;
    real iS;
    real ap;
    real iFp;
    real iSp;
    real d;
    real ff;
    real fs;
    real fcaf;
    real fcas;
    real jca;
    real ffp;
    real fcafp;
    real nca_ss;
    real nca_i;
    real C1;
    real C2;
    real C3;
    real I;
    real O;
    real xs1;
    real xs2;
    real Jrel_np;
    real Jrel_p;

    if (use_adpt_dt) {
        v = sv[0];
        CaMKt = sv[1];
        cass = sv[2];
        nai = sv[3];
        nass = sv[4];
        ki = sv[5];
        kss = sv[6];
        cansr = sv[7];
        cajsr = sv[8];
        cai = sv[9];
        m = sv[10];
        h = sv[11];
        j = sv[12];
        hp = sv[13];
        jp = sv[14];
        mL = sv[15];
        hL = sv[16];
        hLp = sv[17];
        a = sv[18];
        iF = sv[19];
        iS = sv[20];
        ap = sv[21];
        iFp = sv[22];
        iSp = sv[23];
        d = sv[24];
        ff = sv[25];
        fs = sv[26];
        fcaf = sv[27];
        fcas = sv[28];
        jca = sv[29];
        ffp = sv[30];
        fcafp = sv[31];
        nca_ss = sv[32];
        nca_i = sv[33];
        C1 = sv[34];
        C2 = sv[35];
        C3 = sv[36];
        I = sv[37];
        O = sv[38];
        xs1 = sv[39];
        xs2 = sv[40];
        Jrel_np = sv[41];
        Jrel_p = sv[42];
    } else {
        v = *((real *)((char *)sv + pitch * 0) + threadID_);
        CaMKt = *((real *)((char *)sv + pitch * 1) + threadID_);
        cass = *((real *)((char *)sv + pitch * 2) + threadID_);
        nai = *((real *)((char *)sv + pitch * 3) + threadID_);
        nass = *((real *)((char *)sv + pitch * 4) + threadID_);
        ki = *((real *)((char *)sv + pitch * 5) + threadID_);
        kss = *((real *)((char *)sv + pitch * 6) + threadID_);
        cansr = *((real *)((char *)sv + pitch * 7) + threadID_);
        cajsr = *((real *)((char *)sv + pitch * 8) + threadID_);
        cai = *((real *)((char *)sv + pitch * 9) + threadID_);
        m = *((real *)((char *)sv + pitch * 10) + threadID_);
        h = *((real *)((char *)sv + pitch * 11) + threadID_);
        j = *((real *)((char *)sv + pitch * 12) + threadID_);
        hp = *((real *)((char *)sv + pitch * 13) + threadID_);
        jp = *((real *)((char *)sv + pitch * 14) + threadID_);
        mL = *((real *)((char *)sv + pitch * 15) + threadID_);
        hL = *((real *)((char *)sv + pitch * 16) + threadID_);
        hLp = *((real *)((char *)sv + pitch * 17) + threadID_);
        a = *((real *)((char *)sv + pitch * 18) + threadID_);
        iF = *((real *)((char *)sv + pitch * 19) + threadID_);
        iS = *((real *)((char *)sv + pitch * 20) + threadID_);
        ap = *((real *)((char *)sv + pitch * 21) + threadID_);
        iFp = *((real *)((char *)sv + pitch * 22) + threadID_);
        iSp = *((real *)((char *)sv + pitch * 23) + threadID_);
        d = *((real *)((char *)sv + pitch * 24) + threadID_);
        ff = *((real *)((char *)sv + pitch * 25) + threadID_);
        fs = *((real *)((char *)sv + pitch * 26) + threadID_);
        fcaf = *((real *)((char *)sv + pitch * 27) + threadID_);
        fcas = *((real *)((char *)sv + pitch * 28) + threadID_);
        jca = *((real *)((char *)sv + pitch * 29) + threadID_);
        ffp = *((real *)((char *)sv + pitch * 30) + threadID_);
        fcafp = *((real *)((char *)sv + pitch * 31) + threadID_);
        nca_ss = *((real *)((char *)sv + pitch * 32) + threadID_);
        nca_i = *((real *)((char *)sv + pitch * 33) + threadID_);
        C1 = *((real *)((char *)sv + pitch * 34) + threadID_);
        C2 = *((real *)((char *)sv + pitch * 35) + threadID_);
        C3 = *((real *)((char *)sv + pitch * 36) + threadID_);
        I = *((real *)((char *)sv + pitch * 37) + threadID_);
        O = *((real *)((char *)sv + pitch * 38) + threadID_);
        xs1 = *((real *)((char *)sv + pitch * 39) + threadID_);
        xs2 = *((real *)((char *)sv + pitch * 40) + threadID_);
        Jrel_np = *((real *)((char *)sv + pitch * 41) + threadID_);
        Jrel_p = *((real *)((char *)sv + pitch * 42) + threadID_);
    }

    #include "ToRORd_fkatp_mixed_endo_mid_epi_RL.common.c"
}
