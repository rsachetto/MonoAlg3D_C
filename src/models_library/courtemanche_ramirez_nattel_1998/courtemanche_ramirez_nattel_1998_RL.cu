#include "courtemanche_ramirez_nattel_1998.h"
#include <stddef.h>
#include <stdint.h>

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

         *((real * )((char *) sv + pitch * 0) + threadID) = -8.118000e+01f; //V millivolt 
         *((real * )((char *) sv + pitch * 1) + threadID) = 2.908000e-03f; //m dimensionless 
         *((real * )((char *) sv + pitch * 2) + threadID) = 9.649000e-01f; //h dimensionless 
         *((real * )((char *) sv + pitch * 3) + threadID) = 9.775000e-01f; //j dimensionless 
         *((real * )((char *) sv + pitch * 4) + threadID) = 3.043000e-02f; //oa dimensionless 
         *((real * )((char *) sv + pitch * 5) + threadID) = 9.992000e-01f; //oi dimensionless 
         *((real * )((char *) sv + pitch * 6) + threadID) = 4.966000e-03f; //ua dimensionless 
         *((real * )((char *) sv + pitch * 7) + threadID) = 9.986000e-01f; //ui dimensionless 
         *((real * )((char *) sv + pitch * 8) + threadID) = 3.296000e-05f; //xr dimensionless 
         *((real * )((char *) sv + pitch * 9) + threadID) = 1.869000e-02f; //xs dimensionless 
         *((real * )((char *) sv + pitch * 10) + threadID) = 1.367000e-04f; //d dimensionless 
         *((real * )((char *) sv + pitch * 11) + threadID) = 9.996000e-01f; //f dimensionless 
         *((real * )((char *) sv + pitch * 12) + threadID) = 7.755000e-01f; //f_Ca dimensionless 
         *((real * )((char *) sv + pitch * 13) + threadID) = 0.0f; //u dimensionless
         *((real * )((char *) sv + pitch * 14) + threadID) = 1.000000e+00f; //v dimensionless 
         *((real * )((char *) sv + pitch * 15) + threadID) = 9.992000e-01f; //w dimensionless 
         *((real * )((char *) sv + pitch * 16) + threadID) = 1.117000e+01f; //Na_i millimolar 
         *((real * )((char *) sv + pitch * 17) + threadID) = 1.390000e+02f; //K_i millimolar 
         *((real * )((char *) sv + pitch * 18) + threadID) = 1.013000e-04f; //Ca_i millimolar 
         *((real * )((char *) sv + pitch * 19) + threadID) = 1.488000e+00f; //Ca_up millimolar 
         *((real * )((char *) sv + pitch * 20) + threadID) = 1.488000e+00f; //Ca_rel millimolar 

        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * NEQ) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * (NEQ + 1)) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * (NEQ + 2)) + threadID) = 0.0;    // previous dt
        }
	}
}

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int thread_id, real dt, size_t pitch, bool use_adpt_dt) {

    //State variables
    real V_old_; //millivolt
    real m_old_; //dimensionless
    real h_old_; //dimensionless
    real j_old_; //dimensionless
    real oa_old_; //dimensionless
    real oi_old_; //dimensionless
    real ua_old_; //dimensionless
    real ui_old_; //dimensionless
    real xr_old_; //dimensionless
    real xs_old_; //dimensionless
    real d_old_; //dimensionless
    real f_old_; //dimensionless
    real f_Ca_old_; //dimensionless
    real u_old_; //dimensionless
    real v_old_; //dimensionless
    real w_old_; //dimensionless
    real Na_i_old_; //millimolar
    real K_i_old_; //millimolar
    real Ca_i_old_; //millimolar
    real Ca_up_old_; //millimolar
    real Ca_rel_old_; //millimolar


    if(use_adpt_dt) {
        V_old_ =  sv[0];
        m_old_ =  sv[1];
        h_old_ =  sv[2];
        j_old_ =  sv[3];
        oa_old_ =  sv[4];
        oi_old_ =  sv[5];
        ua_old_ =  sv[6];
        ui_old_ =  sv[7];
        xr_old_ =  sv[8];
        xs_old_ =  sv[9];
        d_old_ =  sv[10];
        f_old_ =  sv[11];
        f_Ca_old_ =  sv[12];
        u_old_ =  sv[13];
        v_old_ =  sv[14];
        w_old_ =  sv[15];
        Na_i_old_ =  sv[16];
        K_i_old_ =  sv[17];
        Ca_i_old_ =  sv[18];
        Ca_up_old_ =  sv[19];
        Ca_rel_old_ =  sv[20];
    } else {
        V_old_ =  *((real*)((char*)sv + pitch * 0) + thread_id);
        m_old_ =  *((real*)((char*)sv + pitch * 1) + thread_id);
        h_old_ =  *((real*)((char*)sv + pitch * 2) + thread_id);
        j_old_ =  *((real*)((char*)sv + pitch * 3) + thread_id);
        oa_old_ =  *((real*)((char*)sv + pitch * 4) + thread_id);
        oi_old_ =  *((real*)((char*)sv + pitch * 5) + thread_id);
        ua_old_ =  *((real*)((char*)sv + pitch * 6) + thread_id);
        ui_old_ =  *((real*)((char*)sv + pitch * 7) + thread_id);
        xr_old_ =  *((real*)((char*)sv + pitch * 8) + thread_id);
        xs_old_ =  *((real*)((char*)sv + pitch * 9) + thread_id);
        d_old_ =  *((real*)((char*)sv + pitch * 10) + thread_id);
        f_old_ =  *((real*)((char*)sv + pitch * 11) + thread_id);
        f_Ca_old_ =  *((real*)((char*)sv + pitch * 12) + thread_id);
        u_old_ =  *((real*)((char*)sv + pitch * 13) + thread_id);
        v_old_ =  *((real*)((char*)sv + pitch * 14) + thread_id);
        w_old_ =  *((real*)((char*)sv + pitch * 15) + thread_id);
        Na_i_old_ =  *((real*)((char*)sv + pitch * 16) + thread_id);
        K_i_old_ =  *((real*)((char*)sv + pitch * 17) + thread_id);
        Ca_i_old_ =  *((real*)((char*)sv + pitch * 18) + thread_id);
        Ca_up_old_ =  *((real*)((char*)sv + pitch * 19) + thread_id);
        Ca_rel_old_ =  *((real*)((char*)sv + pitch * 20) + thread_id);
    }

    #include "courtemanche_ramirez_nattel_1998_RL_common.inc.c"

}

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    size_t pitch_h;

    uint8_t use_adpt_dt = (uint8_t)solver->adaptive;

    log_info("Using GPU model implemented in %s\n", __FILE__);

    uint32_t num_volumes = solver->original_num_cells;

    if(use_adpt_dt) {
        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_info("Using Rush-Larsen/Euler model to solve the ODEs\n");
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
			// //computes euler method
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

            for(int n = 0; n < num_steps; ++n) {

                RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, pitch, false);

                // Full Euler
                //for(int i = 0; i < NEQ; i++) {
                //    *((real *)((char *)sv + pitch * i) + sv_id) =
                //        dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
                //}

                // Rush-Larsen/Euler
                *((real *)((char *)sv + pitch * 0) + sv_id) = dt * rDY[0] + *((real *)((char *)sv + pitch * 0) + sv_id);
                *((real *)((char *)sv + pitch * 1) + sv_id) = rDY[1];
                *((real *)((char *)sv + pitch * 2) + sv_id) = rDY[2];
                *((real *)((char *)sv + pitch * 3) + sv_id) = rDY[3];
                *((real *)((char *)sv + pitch * 4) + sv_id) = rDY[4];
                *((real *)((char *)sv + pitch * 5) + sv_id) = rDY[5];
                *((real *)((char *)sv + pitch * 6) + sv_id) = rDY[6];
                *((real *)((char *)sv + pitch * 7) + sv_id) = rDY[7];
                *((real *)((char *)sv + pitch * 8) + sv_id) = rDY[8];
                *((real *)((char *)sv + pitch * 9) + sv_id) = rDY[9];
                *((real *)((char *)sv + pitch * 10) + sv_id) = rDY[10];
                *((real *)((char *)sv + pitch * 11) + sv_id) = rDY[11];
                *((real *)((char *)sv + pitch * 12) + sv_id) = rDY[12];
                *((real *)((char *)sv + pitch * 13) + sv_id) = rDY[13];
                *((real *)((char *)sv + pitch * 14) + sv_id) = rDY[14];
                *((real *)((char *)sv + pitch * 15) + sv_id) = rDY[15];
                *((real *)((char *)sv + pitch * 16) + sv_id) = dt * rDY[16] + *((real *)((char *)sv + pitch * 16) + sv_id);
                *((real *)((char *)sv + pitch * 17) + sv_id) = dt * rDY[17] + *((real *)((char *)sv + pitch * 17) + sv_id);
                *((real *)((char *)sv + pitch * 18) + sv_id) = dt * rDY[18] + *((real *)((char *)sv + pitch * 18) + sv_id);
                *((real *)((char *)sv + pitch * 19) + sv_id) = dt * rDY[19] + *((real *)((char *)sv + pitch * 19) + sv_id);
                *((real *)((char *)sv + pitch * 20) + sv_id) = dt * rDY[20] + *((real *)((char *)sv + pitch * 20) + sv_id);
            }
        } else {
            solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
        }
    }
}

