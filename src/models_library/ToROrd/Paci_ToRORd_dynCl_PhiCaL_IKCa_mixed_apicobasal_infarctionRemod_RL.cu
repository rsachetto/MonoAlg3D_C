#include "Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_RL.h"
#include "../../gpu_utils/gpu_utils.h"
#include <stddef.h>
#include <stdint.h>

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {
    size_t pitch_h;
    uint8_t use_adpt_dt = (uint8_t)solver->adaptive;

    log_info("Using Paci_ToRORd_dynCl_mixed_endo_mid_epi GPU model\n");

    uint32_t num_volumes = solver->original_num_cells;

    if(use_adpt_dt) {
        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_info("Using Euler model to solve the ODEs\n");
    }

    // Execution configuration
    const int GRID = (num_volumes + BLOCK_SIZE - 1) / BLOCK_SIZE;

    size_t size = num_volumes * sizeof(real);

    if(use_adpt_dt)
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ + 3));
    else
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ));

    //Get the extra_data array
    real *extra_data = NULL;
    real *extra_data_device = NULL;

    if(solver->ode_extra_data) {
        extra_data = (real*)solver->ode_extra_data;
        check_cuda_error(cudaMalloc((void **)&extra_data_device, solver->extra_data_size));
        check_cuda_error(cudaMemcpy(extra_data_device, extra_data, solver->extra_data_size, cudaMemcpyHostToDevice));
    }
    kernel_set_model_initial_conditions<<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes, pitch_h, use_adpt_dt, solver->min_dt, extra_data_device);

    check_cuda_error(cudaPeekAtLastError());
    cudaDeviceSynchronize();
    check_cuda_error(cudaFree(extra_data_device));

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

    //Get the extra_data array
    real *extra_data = NULL;
    real *extra_data_device = NULL;
    if(ode_solver->ode_extra_data) {
        extra_data = (real*)ode_solver->ode_extra_data;
        check_cuda_error(cudaMalloc((void **)&extra_data_device, ode_solver->extra_data_size));
        check_cuda_error(cudaMemcpy(extra_data_device, extra_data, ode_solver->extra_data_size, cudaMemcpyHostToDevice));
    } else {
        log_error_and_exit("You need to specify a mask function when using a mixed model!\n");
    }
    real current_scaling[45];
    real *current_scaling_device = NULL;
    for (int i = 1; i < 45; i++) {
        current_scaling[i] = (real) extra_data[5 * num_cells_to_solve + i];
    }
    check_cuda_error(cudaMalloc((void **)&current_scaling_device, sizeof(real) * 45));
    check_cuda_error(cudaMemcpy(current_scaling_device, current_scaling, sizeof(real) * 45, cudaMemcpyHostToDevice));

    solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, ode_solver->pitch, ode_solver->adaptive, ode_solver->abs_tol, ode_solver->rel_tol, ode_solver->max_dt, extra_data_device, current_scaling_device);

    check_cuda_error(cudaPeekAtLastError());

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if(extra_data_device) check_cuda_error(cudaFree(extra_data_device));
}

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt, real *extra_data) {
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int is_paci = (int) extra_data[thread_id];
    int layer = (int) extra_data[thread_id + num_volumes];
    int infarct_zone = (int) extra_data[thread_id + (2*num_volumes)];
    int infarct_stage = (int) extra_data[5 * num_volumes];

    if (thread_id < num_volumes) {
        #include "Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_SS.common.c"
        if (is_paci > 0) {
            *((real * )((char *) sv + pitch * 0) + thread_id)  = v;             // mV;         millivolt
            *((real * )((char *) sv + pitch * 1) + thread_id)  = caSR;          // Ca_SR;     millimolar
            *((real * )((char *) sv + pitch * 2) + thread_id)  = cai;        // Ca_i;      millimolar
            *((real * )((char *) sv + pitch * 3) + thread_id)  = d;       // d;         dimensionless
            *((real * )((char *) sv + pitch * 4) + thread_id)  = f1;          // f1;        dimensionless
            *((real * )((char *) sv + pitch * 5) + thread_id)  = f2;          // f2;        dimensionless
            *((real * )((char *) sv + pitch * 6) + thread_id)  = fCa;          // fCa;       dimensionless
            *((real * )((char *) sv + pitch * 7) + thread_id)  = Xr1;        // Xr1;       dimensionless
            *((real * )((char *) sv + pitch * 8) + thread_id)  = Xr2;          // Xr2;       dimensionless
            *((real * )((char *) sv + pitch * 9) + thread_id)  = Xs;         // Xs;        dimensionless
            *((real * )((char *) sv + pitch *10) + thread_id)  = h;          // h;         dimensionless
            *((real * )((char *) sv + pitch *11) + thread_id)  = j;          // j;         dimensionless
            *((real * )((char *) sv + pitch *12) + thread_id)  = m;         // m;         dimensionless
            *((real * )((char *) sv + pitch *13) + thread_id)  = Xf;          // Xf;        dimensionless
            *((real * )((char *) sv + pitch *14) + thread_id)  = q;          // q;         dimensionless
            *((real * )((char *) sv + pitch *15) + thread_id)  = r;        // r;         dimensionless
            *((real * )((char *) sv + pitch *16) + thread_id)  = nai;           // Na_i;      millimolar
            *((real * )((char *) sv + pitch *17) + thread_id)  = mL;        // m_L;       dimensionless
            *((real * )((char *) sv + pitch *18) + thread_id)  = hL;          // h_L;       dimensionless
            *((real * )((char *) sv + pitch *19) + thread_id)  = RyRa;         // RyRa;      dimensionless
            *((real * )((char *) sv + pitch *20) + thread_id)  = RyRo;       // RyRo;      dimensionless
            *((real * )((char *) sv + pitch *21) + thread_id)  = RyRc;          // RyRc;      dimensionless
            *((real * )((char *) sv + pitch *22) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *23) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *24) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *25) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *26) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *27) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *28) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *29) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *30) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *31) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *32) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *33) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *34) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *35) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *36) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *37) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *38) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *39) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *40) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *41) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *42) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *43) + thread_id)  = 0.0;
            *((real * )((char *) sv + pitch *44) + thread_id)  = 0.0;
        } else if (is_paci == 0) {
            *((real *)((char *)sv + pitch *  0)+thread_id) = v;
            *((real *)((char *)sv + pitch *  1)+thread_id) = nai;
            *((real *)((char *)sv + pitch *  2)+thread_id) = nass;
            *((real *)((char *)sv + pitch *  3)+thread_id) = ki;
            *((real *)((char *)sv + pitch *  4)+thread_id) = kss;
            *((real *)((char *)sv + pitch *  5)+thread_id) = cai;
            *((real *)((char *)sv + pitch *  6)+thread_id) = cass;
            *((real *)((char *)sv + pitch *  7)+thread_id) = cansr;
            *((real *)((char *)sv + pitch *  8)+thread_id) = cajsr;
            *((real *)((char *)sv + pitch *  9)+thread_id) = m;
            *((real *)((char *)sv + pitch * 10)+thread_id) = hp;
            *((real *)((char *)sv + pitch * 11)+thread_id) = h;
            *((real *)((char *)sv + pitch * 12)+thread_id) = j;
            *((real *)((char *)sv + pitch * 13)+thread_id) = jp;
            *((real *)((char *)sv + pitch * 14)+thread_id) = mL;
            *((real *)((char *)sv + pitch * 15)+thread_id) = hL;
            *((real *)((char *)sv + pitch * 16)+thread_id) = hLp;
            *((real *)((char *)sv + pitch * 17)+thread_id) = a;
            *((real *)((char *)sv + pitch * 18)+thread_id) = iF;
            *((real *)((char *)sv + pitch * 19)+thread_id) = iS;
            *((real *)((char *)sv + pitch * 20)+thread_id) = ap;
            *((real *)((char *)sv + pitch * 21)+thread_id) = iFp;
            *((real *)((char *)sv + pitch * 22)+thread_id) = iSp;
            *((real *)((char *)sv + pitch * 23)+thread_id) = d;
            *((real *)((char *)sv + pitch * 24)+thread_id) = ff;
            *((real *)((char *)sv + pitch * 25)+thread_id) = fs;
            *((real *)((char *)sv + pitch * 26)+thread_id) = fcaf;
            *((real *)((char *)sv + pitch * 27)+thread_id) = fcas;
            *((real *)((char *)sv + pitch * 28)+thread_id) = jca;
            *((real *)((char *)sv + pitch * 29)+thread_id) = nca;
            *((real *)((char *)sv + pitch * 30)+thread_id) = nca_i;
            *((real *)((char *)sv + pitch * 31)+thread_id) = ffp;
            *((real *)((char *)sv + pitch * 32)+thread_id) = fcafp;
            *((real *)((char *)sv + pitch * 33)+thread_id) = xs1;
            *((real *)((char *)sv + pitch * 34)+thread_id) = xs2;
            *((real *)((char *)sv + pitch * 35)+thread_id) = Jrel_np;
            *((real *)((char *)sv + pitch * 36)+thread_id) = CaMKt;
            *((real *)((char *)sv + pitch * 37)+thread_id) = ikr_c0;
            *((real *)((char *)sv + pitch * 38)+thread_id) = ikr_c1;
            *((real *)((char *)sv + pitch * 39)+thread_id) = ikr_c2;
            *((real *)((char *)sv + pitch * 40)+thread_id) = ikr_o;
            *((real *)((char *)sv + pitch * 41)+thread_id) = ikr_i;
            *((real *)((char *)sv + pitch * 42)+thread_id) = Jrel_p;
            *((real *)((char *)sv + pitch * 43)+thread_id) = cli;
            *((real *)((char *)sv + pitch * 44)+thread_id) = clss;
        }

        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * (NEQ)) + thread_id) = min_dt;   // dt
            *((real *)((char *)sv + pitch * (NEQ+1)) + thread_id) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * (NEQ+2)) + thread_id) = 0.0;    // previous dt
        }
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents, uint32_t *cells_to_solve,
                            uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt, real *extra_data, real *current_scaling) {

    const real TOLERANCE = 1e-8;
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    int is_paci = (int) extra_data[thread_id];
    int layer = (int) extra_data[thread_id + num_cells_to_solve];
    int infarct_zone = (int) extra_data[thread_id + (2*num_cells_to_solve)];
    real apicobasal = extra_data[thread_id + (3*num_cells_to_solve)];
    int infarct_stage = (int) extra_data[5 * num_cells_to_solve];

    // Each thread solves one cell model
    if(thread_id < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[thread_id];
        else
            sv_id = thread_id;
        if((infarct_zone != 1) || (infarct_stage != 3)) {
            if(!use_adpt) {
                real rDY[NEQ];
                real a[NEQ], b[NEQ];

                bool is_rush_larsen[NEQ];

                //for (int i = 0; i < NEQ; i++)
                //       is_rush_larsen[i] = false;
                if((is_paci > 0) && (is_paci < 4)) {
                  //for (int i = 0; i < NEQ; i++)
                  //     is_rush_larsen[i] = false;
                   for (int i = 0; i < NEQ_PACI; i++)
                        is_rush_larsen[i] = true;
                   for (int i = NEQ_PACI; i < NEQ; i++)
                        is_rush_larsen[i] = false;
                   is_rush_larsen[0] = false;
                   is_rush_larsen[1] = false;
                   is_rush_larsen[2] = false;
                   is_rush_larsen[16] = false;
                } else {
                   for (int i = 0; i < 9; i++)
                       is_rush_larsen[i] = false;
                   for (int i = 9; i < 29; i++)
                       is_rush_larsen[i] = true;
                   is_rush_larsen[29] = false;
                   is_rush_larsen[30] = false;
                   for (int i = 31; i < 36; i++)
                       is_rush_larsen[i] = true;
                   for (int i = 36; i < 42; i++)
                       is_rush_larsen[i] = false;
                   is_rush_larsen[42] = true;
                   is_rush_larsen[43] = false;
                   is_rush_larsen[44] = false;
                }

                for(int n = 0; n < num_steps; ++n) {

                    RHS_gpu(a, b, sv, rDY, stim_currents[thread_id], sv_id, dt, pitch, false, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);

                    for (int i = 0; i < NEQ; i++) {
                        if (is_rush_larsen[i]) {
                            SOLVE_EQUATION_RUSH_LARSEN_GPU(i);
                        }
                        else {
                            SOLVE_EQUATION_EULER_GPU(i);
                        }
                    }
                    // if(is_paci == 1) {
                    //     for(int i = 0; i < NEQ_PACI; i++) {
                    //         *((real *)((char *)sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
                    //     }
                    // } else if(is_paci == 0) {
                    //     for(int i = 0; i < NEQ; i++)
                    //         *((real *)((char *)sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
                    // }
                }
            } else {
                //solve_forward_euler_gpu_adpt(sv, stim_currents[thread_id], cur_time + max_dt, sv_id, pitch, abstol, reltol, dt, max_dt, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
                solve_RL_gpu_adpt(sv, stim_currents[thread_id], cur_time + max_dt, sv_id, pitch, abstol, reltol, dt, max_dt, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
            }
        }
    }
}

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt, int is_paci, int layer, int infarct_zone, int infarct_stage, real apicobasal, real *current_scaling) {
    #define DT *((real *)((char *)sv + pitch * (NEQ)) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * (NEQ+1)) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * (NEQ+2)) + thread_id)

    int num_ODEs;
    if(is_paci > 0) {
        num_ODEs = (int) NEQ_PACI;
    } else {
        num_ODEs = (int) NEQ;
    }

    real rDY[NEQ], a_[NEQ], b_[NEQ], a_new[NEQ], b_new[NEQ];

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

    const real __tiny_ = powf(abstol, 2.0);

    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    for(int i = 0; i < num_ODEs; i++) {
        sv_local[i] = *((real *)((char *)sv + pitch * i) + thread_id);
    }

    RHS_gpu(a_, b_, sv_local, rDY, stim_curr, thread_id, dt, pitch, true, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
    time_new += dt;

    for(int i = 0; i < num_ODEs; i++) {
        _k1__[i] = rDY[i];
    }

  	while(1) {

    		for(int i = 0; i < num_ODEs; i++) {
      			// stores the old variables in a vector
      			edos_old_aux_[i] = sv_local[i];
      			// computes euler method
      			edos_new_euler_[i] = _k1__[i] * dt + edos_old_aux_[i];
      			// steps ahead to compute the rk2 method
      			sv_local[i] = edos_new_euler_[i];
  		  }

  		  time_new += dt;

  		  RHS_gpu(a_new, b_new, sv_local, rDY, stim_curr, thread_id, dt, pitch, true, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
  		  time_new -= dt; // step back

  		  real greatestError = 0.0, auxError = 0.0;

  	  	for(int i = 0; i < num_ODEs; i++) {
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
    		dt = _beta_safety_ * dt * sqrtf(1.0f / greatestError);

    		if(dt < min_dt) {
    		   	dt = min_dt;
    		}	else if(dt > max_dt) {
    			  dt = max_dt;
    		}

    		if(time_new + dt > final_time) {
    			dt = final_time - time_new;
    		}

    		// it doesn't accept the solution or accept and risk a NaN
    		if(greatestError >= 1.0f && dt > min_dt) {
      			// restore the old values to do it again
      			for(int i = 0; i < num_ODEs; i++) {
      				   sv_local[i] = edos_old_aux_[i];
      			}
    		} else {
      			for(int i = 0; i < num_ODEs; i++) {
        				_k_aux__[i] = _k2__[i];
        				_k2__[i] = _k1__[i];
        				_k1__[i] = _k_aux__[i];
      			}

      			for(int i = 0; i < num_ODEs; i++) {
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

    for(int i = 0; i < num_ODEs; i++) {
        *((real *)((char *)sv + pitch * i) + thread_id) = sv_local[i];
    }

    DT = dt;
    TIME_NEW = time_new;
    PREVIOUS_DT = previous_dt;
}

inline __device__ void solve_RL_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt, int is_paci, int layer, int infarct_zone, int infarct_stage, real apicobasal, real *current_scaling) {
    #define DT *((real *)((char *)sv + pitch * (NEQ)) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * (NEQ+1)) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * (NEQ+2)) + thread_id)

    int num_ODEs = NEQ;
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

    const real __tiny_ = powf(abstol, 2.0);

    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    for(int i = 0; i < num_ODEs; i++) {
        sv_local[i] = *((real *)((char *)sv + pitch * i) + thread_id);
    }

    RHS_gpu(a_, b_, sv_local, rDY, stim_curr, thread_id, dt, pitch, true, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
    time_new += dt;

    for(int i = 0; i < num_ODEs; i++) {
        _k1__[i] = rDY[i];
    }

    bool is_rush_larsen[NEQ];
    if((is_paci > 0) && (is_paci < 4)) {
        //for (int i = 0; i < NEQ; i++)
        //     is_rush_larsen[i] = false;
         for (int i = 0; i < NEQ_PACI; i++)
              is_rush_larsen[i] = true;
         for (int i = NEQ_PACI; i < NEQ; i++)
              is_rush_larsen[i] = false;
         is_rush_larsen[0] = false;
         is_rush_larsen[1] = false;
         is_rush_larsen[2] = false;
         is_rush_larsen[16] = false;
    } else {
        for (int i = 0; i < 9; i++)
            is_rush_larsen[i] = false;
        for (int i = 9; i < 29; i++)
            is_rush_larsen[i] = true;
        is_rush_larsen[29] = false;
        is_rush_larsen[30] = false;
        for (int i = 31; i < 36; i++)
            is_rush_larsen[i] = true;
        for (int i = 36; i < 42; i++)
            is_rush_larsen[i] = false;
        is_rush_larsen[42] = true;
        is_rush_larsen[43] = false;
        is_rush_larsen[44] = false;
    }
    // if(is_paci == 1) {
    //     num_ODEs = (int) NEQ_PACI;
    // } else {
    //     num_ODEs = (int) NEQ;
    // }

  	while(1) {

    		for(int i = 0; i < num_ODEs; i++) {
            if (is_rush_larsen[i]) {
                SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(i);
            } else {
                SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(i);
            }
        }

  		  time_new += dt;

  		  RHS_gpu(a_new, b_new, sv_local, rDY, stim_curr, thread_id, dt, pitch, true, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
  		  time_new -= dt; // step back

  		  real greatestError = 0.0, auxError = 0.0;
    		real as, bs, f, y_2nd_order;

  	  	for(int i = 0; i < num_ODEs; i++) {
            if (is_rush_larsen[i]) {
                SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(i);
            } else {
                SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(i);
            }
		    }

    		/// adapt the time step
    		greatestError += __tiny_;
    		previous_dt = dt;

    		/// adapt the time step
    		dt = dt * sqrt(0.5f * reltol / greatestError);

    		if(dt < min_dt) {
    		   	dt = min_dt;
    		}	else if(dt > max_dt) {
    			  dt = max_dt;
    		}

    		if(time_new + dt > final_time) {
    			dt = final_time - time_new;
    		}

    		// it doesn't accept the solution or accept and risk a NaN
    		if(greatestError >= 1.0f && dt > min_dt) {
      			// restore the old values to do it again
      			for(int i = 0; i < num_ODEs; i++) {
      				   sv_local[i] = edos_old_aux_[i];
      			}
    		} else {
      			for(int i = 0; i < num_ODEs; i++) {
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

      			for(int i = 0; i < num_ODEs; i++) {
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

    for(int i = 0; i < num_ODEs; i++) {
        *((real *)((char *)sv + pitch * i) + thread_id) = sv_local[i];
    }

    DT = dt;
    TIME_NEW = time_new;
    PREVIOUS_DT = previous_dt;
}

inline __device__ void RHS_gpu(real *a_, real *b_, real *sv, real *rDY, real stim_current, int thread_id, real dt, size_t pitch, bool use_adpt_dt, int is_paci, int layer, int infarct_zone, int infarct_stage, real apicobasal, real *current_scaling) {
    // State variables
    real v;
    real nai;
    real nass;
    real ki;
    real kss;
    real cai;
    real cass;
    real cansr;
    real cajsr;
    real m;
    real hp;
    real h;
    real j;
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
    real nca;
    real nca_i;
    real ffp;
    real fcafp;
    real xs1;
    real xs2;
    real Jrel_np;
    real CaMKt;
    real ikr_c0;
    real ikr_c1;
    real ikr_c2;
    real ikr_o;
    real ikr_i;
    real Jrel_p;
    real caSR;
    real f1;
    real f2;
    real fCa;
    real Xr1;
    real Xr2;
    real Xs;
    real Xf;
    real q;
    real r;
    real RyRa;
    real RyRo;
    real RyRc;
    real cli;
    real clss;

    real INa_mult = current_scaling[1];
    real IKr_mult = current_scaling[2];
    real ICaL_mult = current_scaling[3];
    real INaL_mult = current_scaling[4];
    real IKs_mult = current_scaling[5];
    real Ito_mult = current_scaling[6];
    real IK1_mult = current_scaling[7];

    real INa_mult_RZ = current_scaling[8];
    real IKr_mult_RZ = current_scaling[9];
    real ICaL_mult_RZ = current_scaling[10];
    real INaL_mult_RZ = current_scaling[11];
    real IKs_mult_RZ = current_scaling[12];
    real Ito_mult_RZ = current_scaling[13];
    real IK1_mult_RZ = current_scaling[14];
    real aCaMK_mult_RZ = current_scaling[15];
    real tau_relp_mult_RZ = current_scaling[16];
    real ICab_mult_RZ = current_scaling[17];
    real Iup_mult_RZ = current_scaling[18];
    real IKCa_mult_RZ = current_scaling[19];
    real IClCa_mult_RZ = current_scaling[20];

    real INa_mult_IZ = current_scaling[21];
    real IKr_mult_IZ = current_scaling[22];
    real ICaL_mult_IZ = current_scaling[23];
    real INaL_mult_IZ = current_scaling[24];
    real IKs_mult_IZ = current_scaling[25];
    real Ito_mult_IZ = current_scaling[26];
    real IK1_mult_IZ = current_scaling[27];
    real aCaMK_mult_IZ = current_scaling[28];
    real tau_relp_mult_IZ = current_scaling[29];
    real ICab_mult_IZ = current_scaling[30];
    real Ko_mult_IZ = current_scaling[31];

    real INa_mult_BZ = current_scaling[32];
    real IKr_mult_BZ = current_scaling[33];
    real ICaL_mult_BZ = current_scaling[34];
    real INaL_mult_BZ = current_scaling[35];
    real IKs_mult_BZ = current_scaling[36];
    real Ito_mult_BZ = current_scaling[37];
    real IK1_mult_BZ = current_scaling[38];
    real aCaMK_mult_BZ = current_scaling[39];
    real tau_relp_mult_BZ = current_scaling[40];
    real ICab_mult_BZ = current_scaling[41];
    real Iup_mult_BZ = current_scaling[42];
    real IKCa_mult_BZ = current_scaling[43];
    real IClCa_mult_BZ = current_scaling[44];

    real atrial_INa = 0.5698;
    real atrial_INaL = 3.4781;
    real atrial_Ito = 3.1498;
    real atrial_IKr = 4.2080;
    real atrial_IKs = 0.6331;
    real atrial_IK1 = 0.3464;
    real atrial_ICaL = 1.7471;
    real atrial_INaCa = 1.5376;
    real atrial_INaK = 1.1698;
    real atrial_Iup = 2.9601;
    real atrial_If = 1.1746;
    real atrial_Irel = 0.0943;
    real atrial_Ileak = 0.1118;

    real nodal_INa = 1.1739;
    real nodal_INaL = 1.8433;
    real nodal_Ito = 2.2774;
    real nodal_IKr = 3.15;
    real nodal_IKs = 0.57;
    real nodal_IK1 = 0.01;
    real nodal_ICaL = 0.95;
    real nodal_INaCa = 1.1768;
    real nodal_INaK = 0.95;
    real nodal_Iup = 0.05;
    real nodal_If = 1.25;
    real nodal_Irel = 1.2028;
    real nodal_Ileak = 2.9065;

    if(is_paci > 0) {
        if(use_adpt_dt) {
            v = sv[0];
            caSR = sv[1];
            cai = sv[2];
            d = sv[3];
            f1 = sv[4];
            f2 = sv[5];
            fCa = sv[6];
            Xr1 = sv[7];
            Xr2 = sv[8];
            Xs = sv[9];
            h = sv[10];
            j = sv[11];
            m = sv[12];
            Xf = sv[13];
            q = sv[14];
            r = sv[15];
            nai = sv[16];
            mL = sv[17];
            hL = sv[18];
            RyRa = sv[19];
            RyRo = sv[20];
            RyRc = sv[21];
        } else {
            v = *((real *)((char *)sv + pitch * 0) + thread_id);        // Membrane variable
            caSR = *((real *)((char *)sv + pitch * 1) + thread_id);
            cai = *((real *)((char *)sv + pitch * 2) + thread_id);      // calcium_dynamics__Ca_i
            d = *((real *)((char *)sv + pitch * 3) + thread_id);        // L type Ca current d gate
            f1 = *((real *)((char *)sv + pitch * 4) + thread_id);       // var_L_type_Ca_current_f1_gate__f1
            f2 = *((real *)((char *)sv + pitch * 5) + thread_id);       // var_L_type_Ca_current_f2_gate__f2
            fCa = *((real *)((char *)sv + pitch * 6) + thread_id);      // L_type_Ca_current__fCa
            Xr1 = *((real *)((char *)sv + pitch * 7) + thread_id);      // Rapid time dependent potassium current Xr1
            Xr2 = *((real *)((char *)sv + pitch * 8) + thread_id);      // Rapid time dependent potassium current Xr2
            Xs = *((real *)((char *)sv + pitch * 9) + thread_id);       // Slow time dependent potassium current Xs
            h = *((real *)((char *)sv + pitch * 10) + thread_id);       // Fast sodium current h gate
            j = *((real *)((char *)sv + pitch * 11) + thread_id);       // Fast sodium current j gate
            m = *((real *)((char *)sv + pitch * 12) + thread_id);       // Fast sodium current m
            Xf = *((real *)((char *)sv + pitch * 13) + thread_id);
            q = *((real *)((char *)sv + pitch * 14) + thread_id);
            r = *((real *)((char *)sv + pitch * 15) + thread_id);       // gating r
            nai = *((real *)((char *)sv + pitch * 16) + thread_id);     // var_sodium_dynamics__Na_i
            mL = *((real *)((char *)sv + pitch * 17) + thread_id);
            hL = *((real *)((char *)sv + pitch * 18) + thread_id);
            RyRa = *((real *)((char *)sv + pitch * 19) + thread_id);
            RyRo = *((real *)((char *)sv + pitch * 20) + thread_id);
            RyRc = *((real *)((char *)sv + pitch * 21) + thread_id);
        }
    } else if (is_paci == 0) {
        if (use_adpt_dt) {
            v = sv[0];
            nai = sv[1];
            nass = sv[2];
            ki = sv[3];
            kss = sv[4];
            cai = sv[5];
            cass = sv[6];
            cansr = sv[7];
            cajsr = sv[8];
            m = sv[9];
            hp = sv[10];
            h = sv[11];
            j = sv[12];
            jp = sv[13];
            mL = sv[14];
            hL = sv[15];
            hLp = sv[16];
            a = sv[17];
            iF = sv[18];
            iS = sv[19];
            ap = sv[20];
            iFp = sv[21];
            iSp = sv[22];
            d = sv[23];
            ff = sv[24];
            fs = sv[25];
            fcaf = sv[26];
            fcas = sv[27];
            jca = sv[28];
            nca = sv[29];
            nca_i = sv[30];
            ffp = sv[31];
            fcafp = sv[32];
            xs1 = sv[33];
            xs2 = sv[34];
            Jrel_np = sv[35];
            CaMKt = sv[36];
            ikr_c0 = sv[37];
            ikr_c1 = sv[38];
            ikr_c2 = sv[39];
            ikr_o = sv[40];
            ikr_i = sv[41];
            Jrel_p = sv[42];
            cli = sv[43];
            clss = sv[44];
        } else {
            v = *((real *)((char *)sv + pitch * 0)+thread_id);
            nai = *((real *)((char *)sv + pitch * 1)+thread_id);
            nass = *((real *)((char *)sv + pitch * 2)+thread_id);
            ki = *((real *)((char *)sv + pitch * 3)+thread_id);
            kss = *((real *)((char *)sv + pitch * 4)+thread_id);
            cai = *((real *)((char *)sv + pitch * 5)+thread_id);
            cass = *((real *)((char *)sv + pitch * 6)+thread_id);
            cansr = *((real *)((char *)sv + pitch * 7)+thread_id);
            cajsr = *((real *)((char *)sv + pitch * 8)+thread_id);
            m = *((real *)((char *)sv + pitch * 9)+thread_id);
            hp = *((real *)((char *)sv + pitch * 10)+thread_id);
            h = *((real *)((char *)sv + pitch * 11)+thread_id);
            j = *((real *)((char *)sv + pitch * 12)+thread_id);
            jp = *((real *)((char *)sv + pitch * 13)+thread_id);
            mL = *((real *)((char *)sv + pitch * 14)+thread_id);
            hL = *((real *)((char *)sv + pitch * 15)+thread_id);
            hLp = *((real *)((char *)sv + pitch * 16)+thread_id);
            a = *((real *)((char *)sv + pitch * 17)+thread_id);
            iF = *((real *)((char *)sv + pitch * 18)+thread_id);
            iS = *((real *)((char *)sv + pitch * 19)+thread_id);
            ap = *((real *)((char *)sv + pitch * 20)+thread_id);
            iFp = *((real *)((char *)sv + pitch * 21)+thread_id);
            iSp = *((real *)((char *)sv + pitch * 22)+thread_id);
            d = *((real *)((char *)sv + pitch * 23)+thread_id);
            ff = *((real *)((char *)sv + pitch * 24)+thread_id);
            fs = *((real *)((char *)sv + pitch * 25)+thread_id);
            fcaf = *((real *)((char *)sv + pitch * 26)+thread_id);
            fcas = *((real *)((char *)sv + pitch * 27)+thread_id);
            jca = *((real *)((char *)sv + pitch * 28)+thread_id);
            nca = *((real *)((char *)sv + pitch * 29)+thread_id);
            nca_i = *((real *)((char *)sv + pitch * 30)+thread_id);
            ffp = *((real *)((char *)sv + pitch * 31)+thread_id);
            fcafp = *((real *)((char *)sv + pitch * 32)+thread_id);
            xs1 = *((real *)((char *)sv + pitch * 33)+thread_id);
            xs2 = *((real *)((char *)sv + pitch * 34)+thread_id);
            Jrel_np = *((real *)((char *)sv + pitch * 35)+thread_id);
            CaMKt = *((real *)((char *)sv + pitch * 36)+thread_id);
            ikr_c0 = *((real *)((char *)sv + pitch * 37)+thread_id);
            ikr_c1 = *((real *)((char *)sv + pitch * 38)+thread_id);
            ikr_c2 = *((real *)((char *)sv + pitch * 39)+thread_id);
            ikr_o = *((real *)((char *)sv + pitch * 40)+thread_id);
            ikr_i = *((real *)((char *)sv + pitch * 41)+thread_id);
            Jrel_p = *((real *)((char *)sv + pitch * 42)+thread_id);
            cli = *((real *)((char *)sv + pitch * 43)+thread_id);
            clss = *((real *)((char *)sv + pitch * 44)+thread_id);
        }
    }

    #include "Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_RL.common.c"
}
