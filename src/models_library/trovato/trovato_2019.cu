#include "trovato_2019.h"
#include <stddef.h>
#include <stdint.h>

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        real STATES[NEQ];
        
        // Steady-state 40 pulses (BCL=1000ms)
        STATES[0] = -86.7099;
        STATES[1] = 0.005431;
        STATES[2] = 0.000104;
        STATES[3] = 8.25533;
        STATES[4] = 8.25502;
        STATES[5] = 8.25503;
        STATES[6] = 143.743;
        STATES[7] = 143.744;
        STATES[8] = 143.744;
        STATES[9] = 4.4e-05;
        STATES[10] = 0.000103;
        STATES[11] = 1.26947;
        STATES[12] = 1.25254;
        STATES[13] = 1.27103;
        STATES[14] = 1.1e-05;
        STATES[15] = 0;
        STATES[16] = 0.006303;
        STATES[17] = 0.789469;
        STATES[18] = 0.789392;
        STATES[19] = 0.791301;
        STATES[20] = 0.580955;
        STATES[21] = 0.791719;
        STATES[22] = 0.000241;
        STATES[23] = 0.463851;
        STATES[24] = 0.239936;
        STATES[25] = 0.000272;
        STATES[26] = 0.646362;
        STATES[27] = 0.98999;
        STATES[28] = 0;
        STATES[29] = 1;
        STATES[30] = 0.926919;
        STATES[31] = 1;
        STATES[32] = 1;
        STATES[33] = 0.999976;
        STATES[34] = 1;
        STATES[35] = 1;
        STATES[36] = 0.005885;
        STATES[37] = 0.000303;
        STATES[38] = 0.994251;
        STATES[39] = 0.000367;
        STATES[40] = 0.566131;
        STATES[41] = 0.189842;
        STATES[42] = 0.000222;
        STATES[43] = 0.233515;
        STATES[44] = 0.997077;
        STATES[45] = 0.471259;

        for (int i = 0; i < NEQ; i++)
            *((real * )((char *) sv + pitch * i) + threadID) = STATES[i];


        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * NEQ) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * (NEQ + 1)) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * (NEQ + 2)) + threadID) = 0.0;    // previous dt
        }
    }
}

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    size_t pitch_h;

    uint8_t use_adpt_dt = (uint8_t)solver->adaptive;

    //TODO: set a model name??
    log_info("Using GPU model implemented in %s\n", __FILE__);

    uint32_t num_volumes = solver->original_num_cells;

    if(use_adpt_dt) {
        log_info("Using Adaptive timestep method to solve the ODEs\n");
    } else {
        log_info("Using Fixed timestep model to solve the ODEs\n");
    }

    // execution configuration
    const int GRID = (num_volumes + BLOCK_SIZE - 1) / BLOCK_SIZE;

    size_t size = num_volumes * sizeof(real);

    if(use_adpt_dt)
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ + 3));
    else
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ));

    kernel_set_model_initial_conditions<<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes, pitch_h, use_adpt_dt, solver->min_dt);

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

    // the array cells to solve is passed when we are using and adaptive mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **)&cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(
            cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve,
                                    num_steps, ode_solver->pitch, ode_solver->adaptive, ode_solver->abs_tol,
                                    ode_solver->rel_tol, ode_solver->max_dt);

    check_cuda_error(cudaPeekAtLastError());

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device)
        check_cuda_error(cudaFree(cells_to_solve_device));
}

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt) {

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

    RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt, pitch, true);
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

		RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt, pitch, true);
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

inline __device__ void solve_forward_euler_gpu_adpt_2(real *sv, real stim_curr, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt) {

    #define DT *((real *)((char *)sv + pitch * (NEQ)) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * (NEQ+1)) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * (NEQ+2)) + thread_id)

    real rDY[NEQ];

    real dt = DT;
    real time_new = TIME_NEW;
    real previous_dt = PREVIOUS_DT;

    real edos_old_aux_[NEQ];
    real edos_new_euler_[NEQ];
    real _k1__[NEQ];
    real _k2__[NEQ];
    real _k_aux__[NEQ];
    real sv_local[NEQ];

    const real TOLERANCE = 1e-8;
    const real _beta_safety_ = 0.8;

    const real __tiny_ = pow(abstol, 2.0);

    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    for(int i = 0; i < NEQ; i++) {
        sv_local[i] = *((real *)((char *)sv + pitch * i) + thread_id);
    }

    RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt, pitch, true);
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

		RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt, pitch, true);
		time_new -= dt; // step back

		real greatestError = 0.0, auxError = 0.0;
		
		for(int i = 0; i < NEQ; i++) {

			// stores the new evaluation
            _k2__[i] = rDY[i];
            real f = (_k1__[i] + _k2__[i]) * 0.5;
            real y_2nd_order = edos_old_aux_[i] + (dt) * f;
            auxError = (fabs(y_2nd_order) < TOLERANCE) ? fabs(edos_new_euler_[i] - TOLERANCE) : fabs( (y_2nd_order - edos_new_euler_[i])/(y_2nd_order) );
			// finds the greatest error between  the steps
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

inline __device__ void solve_rush_larsen_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt) {

    #define DT *((real *)((char *)sv + pitch * (NEQ)) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * (NEQ+1)) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * (NEQ+2)) + thread_id)

    // -------------------------------------------------------------------------------------------
    // MODEL SPECIFIC:
    // set the variables which are non-linear and hodkin-huxley type
    const real TOLERANCE = 1e-8;
    bool is_rush_larsen[NEQ];
    for (int i = 0; i < NEQ; i++) {
        is_rush_larsen[i] = ((i >= 16 && i <= 35) || (i >= 37 && i <= 44)) ? true : false;        
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

    RHS_RL_gpu(a_, b_, sv_local, rDY, stim_curr, thread_id, dt, pitch, true);
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

		RHS_RL_gpu(a_new, b_new, sv_local, rDY, stim_curr, thread_id, dt, pitch, true);
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

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve,
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt,
                          real abstol, real reltol, real max_dt) {
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

            // -------------------------------------------------------------------------------------------
            // MODEL SPECIFIC:
            // set the variables which are non-linear and hodkin-huxley type
            const real TOLERANCE = 1e-8;
            bool is_rush_larsen[NEQ];
            for (int i = 0; i < NEQ; i++) {
                is_rush_larsen[i] = ((i >= 16 && i <= 35) || (i >= 37 && i <= 44)) ? true : false;        
            }
            // -------------------------------------------------------------------------------------------

            for(int n = 0; n < num_steps; ++n) {

                RHS_RL_gpu(a, b, sv, rDY, stim_currents[threadID], sv_id, dt, pitch, false);

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
            //solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
            //solve_forward_euler_gpu_adpt_2(sv, stim_currents[threadID], cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
            solve_rush_larsen_gpu_adpt(sv, stim_currents[threadID], cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
        }
    }
}

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real V;
    real CaMKt;
    real cass;
    real nai;
    real nasl;
    real nass;
    real ki;
    real kss;
    real ksl;
    real cai;
    real casl;
    real cansr;
    real cajsr;
    real cacsr;
    real Jrel1;
    real Jrel2;
    real m;
    real hf;
    real hs;
    real j;
    real hsp;
    real jp;
    real mL;
    real hL;
    real hLp;
    real a;
    real i1;
    real i2;
    real d;
    real ff;
    real fs;
    real fcaf;
    real fcas;
    real jca;
    real ffp;
    real fcafp;
    real nca;
    real b;
    real g;
    real xrf;
    real xrs;
    real xs1;
    real xs2;
    real y;
    real xk1;
    real u;   

    if (use_adpt_dt)
    {
        V = sv[0];
        CaMKt = sv[1];
        cass = sv[2];
        nai = sv[3];
        nasl = sv[4];
        nass = sv[5];
        ki = sv[6];
        kss = sv[7];
        ksl = sv[8];
        cai = sv[9];
        casl = sv[10];
        cansr = sv[11];
        cajsr = sv[12];
        cacsr = sv[13];
        Jrel1 = sv[14];
        Jrel2 = sv[15];
        m = sv[16];
        hf = sv[17];
        hs = sv[18];
        j = sv[19];
        hsp = sv[20];
        jp = sv[21];
        mL = sv[22];
        hL = sv[23];
        hLp = sv[24];
        a = sv[25];
        i1 = sv[26];
        i2 = sv[27];
        d = sv[28];
        ff = sv[29];
        fs = sv[30];
        fcaf = sv[31];
        fcas = sv[32];
        jca = sv[33];
        ffp = sv[34];
        fcafp = sv[35];
        nca = sv[36];
        b = sv[37];
        g = sv[38];
        xrf = sv[39];
        xrs = sv[40];
        xs1 = sv[41];
        xs2 = sv[42];
        y = sv[43];
        xk1 = sv[44];
        u = sv[45];
    }
    else
    {
        V = *((real *)((char *)sv + pitch * 0) + threadID_);
        CaMKt = *((real *)((char *)sv + pitch * 1) + threadID_);
        cass = *((real *)((char *)sv + pitch * 2) + threadID_);
        nai = *((real *)((char *)sv + pitch * 3) + threadID_);
        nasl = *((real *)((char *)sv + pitch * 4) + threadID_);
        nass = *((real *)((char *)sv + pitch * 5) + threadID_);
        ki = *((real *)((char *)sv + pitch * 6) + threadID_);
        kss = *((real *)((char *)sv + pitch * 7) + threadID_);
        ksl = *((real *)((char *)sv + pitch * 8) + threadID_);
        cai = *((real *)((char *)sv + pitch * 9) + threadID_);
        casl = *((real *)((char *)sv + pitch * 10) + threadID_);
        cansr = *((real *)((char *)sv + pitch * 11) + threadID_);
        cajsr = *((real *)((char *)sv + pitch * 12) + threadID_);
        cacsr = *((real *)((char *)sv + pitch * 13) + threadID_);
        Jrel1 = *((real *)((char *)sv + pitch * 14) + threadID_);
        Jrel2 = *((real *)((char *)sv + pitch * 15) + threadID_);
        m = *((real *)((char *)sv + pitch * 16) + threadID_);
        hf = *((real *)((char *)sv + pitch * 17) + threadID_);
        hs = *((real *)((char *)sv + pitch * 18) + threadID_);
        j = *((real *)((char *)sv + pitch * 19) + threadID_);
        hsp = *((real *)((char *)sv + pitch * 20) + threadID_);
        jp = *((real *)((char *)sv + pitch * 21) + threadID_);
        mL = *((real *)((char *)sv + pitch * 22) + threadID_);
        hL = *((real *)((char *)sv + pitch * 23) + threadID_);
        hLp = *((real *)((char *)sv + pitch * 24) + threadID_);
        a = *((real *)((char *)sv + pitch * 25) + threadID_);
        i1 = *((real *)((char *)sv + pitch * 26) + threadID_);
        i2 = *((real *)((char *)sv + pitch * 27) + threadID_);
        d = *((real *)((char *)sv + pitch * 28) + threadID_);
        ff = *((real *)((char *)sv + pitch * 29) + threadID_);
        fs = *((real *)((char *)sv + pitch * 30) + threadID_);
        fcaf = *((real *)((char *)sv + pitch * 31) + threadID_);
        fcas = *((real *)((char *)sv + pitch * 32) + threadID_);
        jca = *((real *)((char *)sv + pitch * 33) + threadID_);
        ffp = *((real *)((char *)sv + pitch * 34) + threadID_);
        fcafp = *((real *)((char *)sv + pitch * 35) + threadID_);
        nca = *((real *)((char *)sv + pitch * 36) + threadID_);
        b = *((real *)((char *)sv + pitch * 37) + threadID_);
        g = *((real *)((char *)sv + pitch * 38) + threadID_);
        xrf = *((real *)((char *)sv + pitch * 39) + threadID_);
        xrs = *((real *)((char *)sv + pitch * 40) + threadID_);
        xs1 = *((real *)((char *)sv + pitch * 41) + threadID_);
        xs2 = *((real *)((char *)sv + pitch * 42) + threadID_);
        y = *((real *)((char *)sv + pitch * 43) + threadID_);
        xk1 = *((real *)((char *)sv + pitch * 44) + threadID_);
        u = *((real *)((char *)sv + pitch * 45) + threadID_);
    }

    #include "trovato_2019_common.inc.c"
}

inline __device__ void RHS_RL_gpu(real *a_, real *b_, real *sv, real *rDY_, real stim_current, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real V;
    real CaMKt;
    real cass;
    real nai;
    real nasl;
    real nass;
    real ki;
    real kss;
    real ksl;
    real cai;
    real casl;
    real cansr;
    real cajsr;
    real cacsr;
    real Jrel1;
    real Jrel2;
    real m;
    real hf;
    real hs;
    real j;
    real hsp;
    real jp;
    real mL;
    real hL;
    real hLp;
    real a;
    real i1;
    real i2;
    real d;
    real ff;
    real fs;
    real fcaf;
    real fcas;
    real jca;
    real ffp;
    real fcafp;
    real nca;
    real b;
    real g;
    real xrf;
    real xrs;
    real xs1;
    real xs2;
    real y;
    real xk1;
    real u;   

    if (use_adpt_dt)
    {
        V = sv[0];
        CaMKt = sv[1];
        cass = sv[2];
        nai = sv[3];
        nasl = sv[4];
        nass = sv[5];
        ki = sv[6];
        kss = sv[7];
        ksl = sv[8];
        cai = sv[9];
        casl = sv[10];
        cansr = sv[11];
        cajsr = sv[12];
        cacsr = sv[13];
        Jrel1 = sv[14];
        Jrel2 = sv[15];
        m = sv[16];
        hf = sv[17];
        hs = sv[18];
        j = sv[19];
        hsp = sv[20];
        jp = sv[21];
        mL = sv[22];
        hL = sv[23];
        hLp = sv[24];
        a = sv[25];
        i1 = sv[26];
        i2 = sv[27];
        d = sv[28];
        ff = sv[29];
        fs = sv[30];
        fcaf = sv[31];
        fcas = sv[32];
        jca = sv[33];
        ffp = sv[34];
        fcafp = sv[35];
        nca = sv[36];
        b = sv[37];
        g = sv[38];
        xrf = sv[39];
        xrs = sv[40];
        xs1 = sv[41];
        xs2 = sv[42];
        y = sv[43];
        xk1 = sv[44];
        u = sv[45];
    }
    else
    {
        V = *((real *)((char *)sv + pitch * 0) + threadID_);
        CaMKt = *((real *)((char *)sv + pitch * 1) + threadID_);
        cass = *((real *)((char *)sv + pitch * 2) + threadID_);
        nai = *((real *)((char *)sv + pitch * 3) + threadID_);
        nasl = *((real *)((char *)sv + pitch * 4) + threadID_);
        nass = *((real *)((char *)sv + pitch * 5) + threadID_);
        ki = *((real *)((char *)sv + pitch * 6) + threadID_);
        kss = *((real *)((char *)sv + pitch * 7) + threadID_);
        ksl = *((real *)((char *)sv + pitch * 8) + threadID_);
        cai = *((real *)((char *)sv + pitch * 9) + threadID_);
        casl = *((real *)((char *)sv + pitch * 10) + threadID_);
        cansr = *((real *)((char *)sv + pitch * 11) + threadID_);
        cajsr = *((real *)((char *)sv + pitch * 12) + threadID_);
        cacsr = *((real *)((char *)sv + pitch * 13) + threadID_);
        Jrel1 = *((real *)((char *)sv + pitch * 14) + threadID_);
        Jrel2 = *((real *)((char *)sv + pitch * 15) + threadID_);
        m = *((real *)((char *)sv + pitch * 16) + threadID_);
        hf = *((real *)((char *)sv + pitch * 17) + threadID_);
        hs = *((real *)((char *)sv + pitch * 18) + threadID_);
        j = *((real *)((char *)sv + pitch * 19) + threadID_);
        hsp = *((real *)((char *)sv + pitch * 20) + threadID_);
        jp = *((real *)((char *)sv + pitch * 21) + threadID_);
        mL = *((real *)((char *)sv + pitch * 22) + threadID_);
        hL = *((real *)((char *)sv + pitch * 23) + threadID_);
        hLp = *((real *)((char *)sv + pitch * 24) + threadID_);
        a = *((real *)((char *)sv + pitch * 25) + threadID_);
        i1 = *((real *)((char *)sv + pitch * 26) + threadID_);
        i2 = *((real *)((char *)sv + pitch * 27) + threadID_);
        d = *((real *)((char *)sv + pitch * 28) + threadID_);
        ff = *((real *)((char *)sv + pitch * 29) + threadID_);
        fs = *((real *)((char *)sv + pitch * 30) + threadID_);
        fcaf = *((real *)((char *)sv + pitch * 31) + threadID_);
        fcas = *((real *)((char *)sv + pitch * 32) + threadID_);
        jca = *((real *)((char *)sv + pitch * 33) + threadID_);
        ffp = *((real *)((char *)sv + pitch * 34) + threadID_);
        fcafp = *((real *)((char *)sv + pitch * 35) + threadID_);
        nca = *((real *)((char *)sv + pitch * 36) + threadID_);
        b = *((real *)((char *)sv + pitch * 37) + threadID_);
        g = *((real *)((char *)sv + pitch * 38) + threadID_);
        xrf = *((real *)((char *)sv + pitch * 39) + threadID_);
        xrs = *((real *)((char *)sv + pitch * 40) + threadID_);
        xs1 = *((real *)((char *)sv + pitch * 41) + threadID_);
        xs2 = *((real *)((char *)sv + pitch * 42) + threadID_);
        y = *((real *)((char *)sv + pitch * 43) + threadID_);
        xk1 = *((real *)((char *)sv + pitch * 44) + threadID_);
        u = *((real *)((char *)sv + pitch * 45) + threadID_);
    }

    #include "trovato_2019_RL_common.inc.c"
}
