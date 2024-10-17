#include "ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.h"
#include <stddef.h>
#include <stdint.h>

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {
        for (int i = 0; i < NEQ; i++) {
            // Default initial condition from Matlab (endocardium cell)
            *((real * )((char *) sv + pitch * 0) + threadID) = -88.6369922306458;
            *((real * )((char *) sv + pitch * 1) + threadID) = 11.8973412949238;
            *((real * )((char *) sv + pitch * 2) + threadID) = 11.897661047085;
            *((real * )((char *) sv + pitch * 3) + threadID) = 141.234464714982;
            *((real * )((char *) sv + pitch * 4) + threadID) = 141.234423402713;
            *((real * )((char *) sv + pitch * 5) + threadID) = 7.26747296460659e-05;
            *((real * )((char *) sv + pitch * 6) + threadID) = 6.33786975780735e-05;
            *((real * )((char *) sv + pitch * 7) + threadID) = 1.5326530637197;
            *((real * )((char *) sv + pitch * 8) + threadID) = 1.53394579180493;
            *((real * )((char *) sv + pitch * 9) + threadID) = 0.000828007761976018;
            *((real * )((char *) sv + pitch * 10) + threadID) = 0.666527193684116;
            *((real * )((char *) sv + pitch * 11) + threadID) = 0.826020806005678;
            *((real * )((char *) sv + pitch * 12) + threadID) = 0.826055985895856;
            *((real * )((char *) sv + pitch * 13) + threadID) = 0.825850881115628;
            *((real * )((char *) sv + pitch * 14) + threadID) = 0.000166868626513013;
            *((real * )((char *) sv + pitch * 15) + threadID) = 0.522830604669169;
            *((real * )((char *) sv + pitch * 16) + threadID) = 0.285969584294187;
            *((real * )((char *) sv + pitch * 17) + threadID) = 0.000959137028030184;
            *((real * )((char *) sv + pitch * 18) + threadID) = 0.999601150012565;
            *((real * )((char *) sv + pitch * 19) + threadID) = 0.5934016398361;
            *((real * )((char *) sv + pitch * 20) + threadID) = 0.000488696137242056;
            *((real * )((char *) sv + pitch * 21) + threadID) = 0.999601147267179;
            *((real * )((char *) sv + pitch * 22) + threadID) = 0.654668660159696;
            *((real * )((char *) sv + pitch * 23) + threadID) = 9.50007519781516e-32;
            *((real * )((char *) sv + pitch * 24) + threadID) = 0.999999992317577;
            *((real * )((char *) sv + pitch * 25) + threadID) = 0.939258048397962;
            *((real * )((char *) sv + pitch * 26) + threadID) = 0.999999992317557;
            *((real * )((char *) sv + pitch * 27) + threadID) = 0.999898379647465;
            *((real * )((char *) sv + pitch * 28) + threadID) = 0.99997825156004;
            *((real * )((char *) sv + pitch * 29) + threadID) = 0.000444816183420527;
            *((real * )((char *) sv + pitch * 30) + threadID) = 0.000755072490632667;
            *((real * )((char *) sv + pitch * 31) + threadID) = 0.999999992318446;
            *((real * )((char *) sv + pitch * 32) + threadID) = 0.999999992318445;
            *((real * )((char *) sv + pitch * 33) + threadID) = 0.24240468344952;
            *((real * )((char *) sv + pitch * 34) + threadID) = 0.000179537726989804;
            *((real * )((char *) sv + pitch * 35) + threadID) = -6.88308558109975e-25;
            *((real * )((char *) sv + pitch * 36) + threadID) = 0.0111749845355653;
            *((real * )((char *) sv + pitch * 37) + threadID) = 0.998036620213316;
            *((real * )((char *) sv + pitch * 38) + threadID) = 0.000858801779013532;
            *((real * )((char *) sv + pitch * 39) + threadID) = 0.000709744678350176;
            *((real * )((char *) sv + pitch * 40) + threadID) = 0.000381261722195702;
            *((real * )((char *) sv + pitch * 41) + threadID) = 1.35711566929992e-05;
            *((real * )((char *) sv + pitch * 42) + threadID) = 2.30252452954649e-23;
        }
            
        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * 43) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 44) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 45) + threadID) = 0.0;    // previous dt
        }
    }
}

__global__ void kernel_set_model_initial_conditions_endo_mid_epi(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt,\
                                                real *initial_endo, real *initial_epi, real *initial_mid, real *transmurality, real *sf_Iks) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {
        for (int i = 0; i < NEQ; i++) {
            if (transmurality[threadID] == ENDO)
                *((real * )((char *) sv + pitch * i) + threadID) = initial_endo[i];
            else if (transmurality[threadID] == EPI)
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
    real *transmurality = NULL;
	real *sf_Iks = NULL;
    real *initial_conditions_endo_device = NULL;
    real *initial_conditions_epi_device = NULL;
    real *initial_conditions_mid_device = NULL;
    real *transmurality_device = NULL;
	real *sf_Iks_device = NULL;

    if(solver->ode_extra_data) {
        struct extra_data_for_torord_gksgkrtjca_twave *extra_data = (struct extra_data_for_torord_gksgkrtjca_twave*)solver->ode_extra_data;
        initial_conditions_endo = extra_data->initial_ss_endo;
        initial_conditions_epi = extra_data->initial_ss_epi;
        initial_conditions_mid = extra_data->initial_ss_mid;
        transmurality = extra_data->transmurality;
		sf_Iks = extra_data->sf_IKs;
        check_cuda_error(cudaMalloc((void **)&initial_conditions_endo_device, sizeof(real)*NEQ));
        check_cuda_error(cudaMemcpy(initial_conditions_endo_device, initial_conditions_endo, sizeof(real)*NEQ, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&initial_conditions_epi_device, sizeof(real)*NEQ));
        check_cuda_error(cudaMemcpy(initial_conditions_epi_device, initial_conditions_epi, sizeof(real)*NEQ, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&initial_conditions_mid_device, sizeof(real)*NEQ));
        check_cuda_error(cudaMemcpy(initial_conditions_mid_device, initial_conditions_mid, sizeof(real)*NEQ, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&transmurality_device, sizeof(real)*num_volumes));
        check_cuda_error(cudaMemcpy(transmurality_device, transmurality, sizeof(real)*num_volumes, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&sf_Iks_device, sizeof(real)*num_volumes));
        check_cuda_error(cudaMemcpy(sf_Iks_device, sf_Iks, sizeof(real)*num_volumes, cudaMemcpyHostToDevice));
    }
    else {
        log_info("[INFO] You should supply a mask function to tag the cells when using this mixed model!\n");
        log_info("[INFO] Considering all cells ENDO!\n");
    }

    if (solver->ode_extra_data) {
        kernel_set_model_initial_conditions_endo_mid_epi<<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes, pitch_h, use_adpt_dt, solver->min_dt,\
                                                            initial_conditions_endo_device, initial_conditions_epi_device, initial_conditions_mid_device,\
                                                            transmurality_device, sf_Iks_device);
    }
    else {
        kernel_set_model_initial_conditions<<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes, pitch_h, use_adpt_dt, solver->min_dt);
    }
    

    check_cuda_error(cudaPeekAtLastError());
    cudaDeviceSynchronize();

    check_cuda_error(cudaFree(initial_conditions_endo_device));
    check_cuda_error(cudaFree(initial_conditions_epi_device));
    check_cuda_error(cudaFree(initial_conditions_mid_device));
    check_cuda_error(cudaFree(transmurality_device));

    return pitch_h;
}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_gpu) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;
    bool has_extra_params = (ode_solver->ode_extra_data) ? true : false;

    // execution configuration
    const int GRID = ((int)num_cells_to_solve + BLOCK_SIZE - 1) / BLOCK_SIZE;

    size_t stim_currents_size = sizeof(real) * num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t) * num_cells_to_solve;

    real *stims_currents_device = NULL;
    check_cuda_error(cudaMalloc((void **)&stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    // the array cells to solve is passed when we are using and adaptive mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **)&cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    uint32_t num_volumes = ode_solver->original_num_cells;
    real *transmurality = NULL;
    real *transmurality_device = NULL;
	real *sf_Iks = NULL;
	real *sf_Iks_device = NULL;
    int num_extra_parameters = 20;
    real extra_par[num_extra_parameters];
    real *extra_par_device = NULL;
    
    // Get the extra data array if exists
    // Transmurality and sf_IKs mapping defined on 'extra_data' function
    if(ode_solver->ode_extra_data) {
        struct extra_data_for_torord_gksgkrtjca_twave *extra_data = (struct extra_data_for_torord_gksgkrtjca_twave*)ode_solver->ode_extra_data;
        extra_par[0]  = extra_data->INa_Multiplier; 
        extra_par[1]  = extra_data->INaL_Multiplier;
        extra_par[2]  = extra_data->INaCa_Multiplier;
        extra_par[3]  = extra_data->INaK_Multiplier;
        extra_par[4]  = extra_data->INab_Multiplier; 
        extra_par[5]  = extra_data->Ito_Multiplier;
        extra_par[6]  = extra_data->IKr_Multiplier; 
        extra_par[7]  = extra_data->IKs_Multiplier; 
        extra_par[8]  = extra_data->IK1_Multiplier;
        extra_par[9]  = extra_data->IKb_Multiplier;
        extra_par[10]  = extra_data->IKCa_Multiplier;
        extra_par[11] = extra_data->ICaL_Multiplier;  
        extra_par[12] = extra_data->ICab_Multiplier;  
        extra_par[13] = extra_data->IpCa_Multiplier;
        extra_par[14] = extra_data->ICaCl_Multiplier; 
        extra_par[15] = extra_data->IClb_Multiplier;
        extra_par[16] = extra_data->Jrel_Multiplier;
        extra_par[17] = extra_data->Jup_Multiplier;
        extra_par[18] = extra_data->aCaMK_Multiplier;
        extra_par[19] = extra_data->taurelp_Multiplier;
        sf_Iks = extra_data->sf_IKs;
        transmurality = extra_data->transmurality;

        check_cuda_error(cudaMalloc((void **)&transmurality_device, sizeof(real)*num_volumes));
        check_cuda_error(cudaMemcpy(transmurality_device, transmurality, sizeof(real)*num_volumes, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&sf_Iks_device, sizeof(real)*num_volumes));
        check_cuda_error(cudaMemcpy(sf_Iks_device, sf_Iks, sizeof(real)*num_volumes, cudaMemcpyHostToDevice));
        check_cuda_error(cudaMalloc((void **)&extra_par_device, sizeof(real)*num_extra_parameters));
        check_cuda_error(cudaMemcpy(extra_par_device, extra_par, sizeof(real)*num_extra_parameters, cudaMemcpyHostToDevice));
    }
    // No [extra_data] section, we consider all cells ENDO!
    else {

        // Default: initialize all current modifiers
        for (uint32_t i = 0; i < num_extra_parameters; i++) {
            if (i == 10)
                extra_par[i] = 0.0;
            else 
                extra_par[i] = 1.0;
        }
        check_cuda_error(cudaMalloc((void **)&extra_par_device, sizeof(real)*num_extra_parameters));
        check_cuda_error(cudaMemcpy(extra_par_device, extra_par, sizeof(real)*num_extra_parameters, cudaMemcpyHostToDevice));
    }

    // Call the kernel function to solve the cellular model on the GPU
    solve_endo_mid_epi_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, transmurality_device, sf_Iks_device, extra_par_device,\
                                    num_cells_to_solve, num_steps, ode_solver->pitch, ode_solver->adaptive, ode_solver->abs_tol, ode_solver->rel_tol, ode_solver->max_dt, has_extra_params);

    check_cuda_error(cudaPeekAtLastError());

    if (stims_currents_device) check_cuda_error(cudaFree(stims_currents_device));
    if (cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if (transmurality_device) check_cuda_error(cudaFree(transmurality_device));
	if (sf_Iks_device) check_cuda_error(cudaFree(sf_Iks_device));
    if (extra_par_device) check_cuda_error(cudaFree(extra_par_device));
}

__global__ void solve_endo_mid_epi_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, real *transmurality, real *sf_Iks, real *extra_params,\
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt, bool has_extra_params) {
    const real TOLERANCE = 1e-8;
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

            for(int n = 0; n < num_steps; ++n) {

                // [!] Transmurality and sf_IKs are given by the columns of ALG mesh file 
                if (has_extra_params) {
                    RHS_RL_gpu(a, b, sv, rDY, stim_currents[threadID], transmurality[threadID], sf_Iks[threadID], extra_params, sv_id, dt, pitch, false);
                }
                // Default: all cells are ENDO and sf_IKs equals to 1
                else {
                    RHS_RL_gpu(a, b, sv, rDY, stim_currents[threadID], 0.0, 1.0, extra_params, sv_id, dt, pitch, false);
                }                

                // Solve variables based on its type:
                //  Non-linear = Euler
                //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
                SOLVE_EQUATION_EULER_GPU(0);        // v        
                SOLVE_EQUATION_EULER_GPU(1);        // nai
                SOLVE_EQUATION_EULER_GPU(2);        // nass
                SOLVE_EQUATION_EULER_GPU(3);        // ki   
                SOLVE_EQUATION_EULER_GPU(4);        // kss
                SOLVE_EQUATION_EULER_GPU(5);        // cai
                SOLVE_EQUATION_EULER_GPU(6);        // cass 
                SOLVE_EQUATION_EULER_GPU(7);        // cansr
                SOLVE_EQUATION_EULER_GPU(8);        // cajsr
                SOLVE_EQUATION_RUSH_LARSEN_GPU(9);  // m
                SOLVE_EQUATION_RUSH_LARSEN_GPU(10); // hp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(11); // h
                SOLVE_EQUATION_RUSH_LARSEN_GPU(12); // j
                SOLVE_EQUATION_RUSH_LARSEN_GPU(13); // jp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(14); // mL
                SOLVE_EQUATION_RUSH_LARSEN_GPU(15); // hL
                SOLVE_EQUATION_RUSH_LARSEN_GPU(16); // hLp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(17); // a
                SOLVE_EQUATION_RUSH_LARSEN_GPU(18); // iF
                SOLVE_EQUATION_RUSH_LARSEN_GPU(19); // iS
                SOLVE_EQUATION_RUSH_LARSEN_GPU(20); // ap
                SOLVE_EQUATION_RUSH_LARSEN_GPU(21); // iFp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(22); // iSp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(23); // d
                SOLVE_EQUATION_RUSH_LARSEN_GPU(24); // ff
                SOLVE_EQUATION_RUSH_LARSEN_GPU(25); // fs
                SOLVE_EQUATION_RUSH_LARSEN_GPU(26); // fcaf
                SOLVE_EQUATION_RUSH_LARSEN_GPU(27); // fcas
                SOLVE_EQUATION_RUSH_LARSEN_GPU(28); // jca
                SOLVE_EQUATION_EULER_GPU(29);       // nca
                SOLVE_EQUATION_EULER_GPU(30);       // nca_i
                SOLVE_EQUATION_RUSH_LARSEN_GPU(31); // ffp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(32); // fcafp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(33); // xs1
                SOLVE_EQUATION_RUSH_LARSEN_GPU(34); // xs2
                SOLVE_EQUATION_RUSH_LARSEN_GPU(35); // Jrel_np
                SOLVE_EQUATION_EULER_GPU(36);        // CaMKt 
                SOLVE_EQUATION_EULER_GPU(37);       // ikr_c0
                SOLVE_EQUATION_EULER_GPU(38);       // ikr_c1
                SOLVE_EQUATION_EULER_GPU(39);       // ikr_c2
                SOLVE_EQUATION_EULER_GPU(40);       // ikr_o
                SOLVE_EQUATION_EULER_GPU(41);       // ikr_i
                SOLVE_EQUATION_RUSH_LARSEN_GPU(42); // Jrel_p
            }
        } else {
            solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], transmurality[threadID], sf_Iks[threadID], extra_params, cur_time + max_dt, sv_id, pitch, abstol, reltol, dt, max_dt);
        }
    }
}

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real transmurality, real sf_Iks, real *extra_params, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt) {

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

    RHS_gpu(sv_local, rDY, stim_curr, transmurality, sf_Iks, extra_params, thread_id, dt, pitch, true);
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

		RHS_gpu(sv_local, rDY, stim_curr, transmurality, sf_Iks, extra_params, thread_id, dt, pitch, true);
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

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, real transmurality, real sf_Iks, real *extra_params, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {
    
    // Current modifiers
    real INa_Multiplier = extra_params[0];   
    real INaL_Multiplier = extra_params[1];  
    real INaCa_Multiplier = extra_params[2];  
    real INaK_Multiplier = extra_params[3];  
    real INab_Multiplier = extra_params[4];   
    real Ito_Multiplier = extra_params[5];  
    real IKr_Multiplier = extra_params[6];   
    real IKs_Multiplier = extra_params[7];   
    real IK1_Multiplier = extra_params[8];  
    real IKb_Multiplier = extra_params[9];  
    real IKCa_Multiplier = extra_params[10];  
    real ICaL_Multiplier = extra_params[11];   
    real ICab_Multiplier = extra_params[12];   
    real IpCa_Multiplier = extra_params[13]; 
    real ICaCl_Multiplier = extra_params[14];  
    real IClb_Multiplier = extra_params[15]; 
    real Jrel_Multiplier = extra_params[16]; 
    real Jup_Multiplier = extra_params[17]; 
    real aCaMK_Multiplier = extra_params[18]; 
    real taurelp_Multiplier = extra_params[19];

    // Get the celltype for the current cell
    real celltype = transmurality;
    
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
    real nca;
    real nca_i;
    real ikr_c0;
    real ikr_c1;
    real ikr_c2;
    real ikr_i;
    real ikr_o;
    real xs1;
    real xs2;
    real Jrel_np;
    real Jrel_p;

    if (use_adpt_dt) {
        v       = sv[0];
        nai     = sv[1];
        nass    = sv[2];
        ki      = sv[3];
        kss     = sv[4];
        cai     = sv[5];
        cass    = sv[6];
        cansr   = sv[7];
        cajsr   = sv[8];
        m       = sv[9];
        hp      = sv[10];
        h       = sv[11];
        j       = sv[12];
        jp      = sv[13];
        mL      = sv[14];
        hL      = sv[15];
        hLp     = sv[16];
        a       = sv[17];
        iF      = sv[18];
        iS      = sv[19];
        ap      = sv[20];
        iFp     = sv[21];
        iSp     = sv[22];
        d       = sv[23];
        ff      = sv[24];
        fs      = sv[25];
        fcaf    = sv[26];
        fcas    = sv[27];
        jca     = sv[28];
        nca     = sv[29];
        nca_i   = sv[30];
        ffp     = sv[31];
        fcafp   = sv[32];
        xs1     = sv[33];
        xs2     = sv[34];
        Jrel_np = sv[35];
        CaMKt   = sv[36];
        ikr_c0  = sv[37];
        ikr_c1  = sv[38];
        ikr_c2  = sv[39];
        ikr_o   = sv[40];
        ikr_i   = sv[41];
        Jrel_p  = sv[42];
    } else {
        v       = *((real *)((char *)sv + pitch * 0) + threadID_);
        nai     = *((real *)((char *)sv + pitch * 1) + threadID_);
        nass    = *((real *)((char *)sv + pitch * 2) + threadID_);
        ki      = *((real *)((char *)sv + pitch * 3) + threadID_);
        kss     = *((real *)((char *)sv + pitch * 4) + threadID_);
        cai     = *((real *)((char *)sv + pitch * 5) + threadID_);
        cass    = *((real *)((char *)sv + pitch * 6) + threadID_);
        cansr   = *((real *)((char *)sv + pitch * 7) + threadID_);
        cajsr   = *((real *)((char *)sv + pitch * 8) + threadID_);
        m       = *((real *)((char *)sv + pitch * 9) + threadID_);
        hp      = *((real *)((char *)sv + pitch * 10) + threadID_);
        h       = *((real *)((char *)sv + pitch * 11) + threadID_);
        j       = *((real *)((char *)sv + pitch * 12) + threadID_);
        jp      = *((real *)((char *)sv + pitch * 13) + threadID_);
        mL      = *((real *)((char *)sv + pitch * 14) + threadID_);
        hL      = *((real *)((char *)sv + pitch * 15) + threadID_);
        hLp     = *((real *)((char *)sv + pitch * 16) + threadID_);
        a       = *((real *)((char *)sv + pitch * 17) + threadID_);
        iF      = *((real *)((char *)sv + pitch * 18) + threadID_);
        iS      = *((real *)((char *)sv + pitch * 19) + threadID_);
        ap      = *((real *)((char *)sv + pitch * 20) + threadID_);
        iFp     = *((real *)((char *)sv + pitch * 21) + threadID_);
        iSp     = *((real *)((char *)sv + pitch * 22) + threadID_);
        d       = *((real *)((char *)sv + pitch * 23) + threadID_);
        ff      = *((real *)((char *)sv + pitch * 24) + threadID_);
        fs      = *((real *)((char *)sv + pitch * 25) + threadID_);
        fcaf    = *((real *)((char *)sv + pitch * 26) + threadID_);
        fcas    = *((real *)((char *)sv + pitch * 27) + threadID_);
        jca     = *((real *)((char *)sv + pitch * 28) + threadID_);
        nca     = *((real *)((char *)sv + pitch * 29) + threadID_);
        nca_i   = *((real *)((char *)sv + pitch * 30) + threadID_);
        ffp     = *((real *)((char *)sv + pitch * 31) + threadID_);
        fcafp   = *((real *)((char *)sv + pitch * 32) + threadID_);
        xs1     = *((real *)((char *)sv + pitch * 33) + threadID_);
        xs2     = *((real *)((char *)sv + pitch * 34) + threadID_);
        Jrel_np = *((real *)((char *)sv + pitch * 35) + threadID_);
        CaMKt   = *((real *)((char *)sv + pitch * 36) + threadID_);
        ikr_c0  = *((real *)((char *)sv + pitch * 37) + threadID_);
        ikr_c1  = *((real *)((char *)sv + pitch * 38) + threadID_);
        ikr_c2  = *((real *)((char *)sv + pitch * 39) + threadID_);
        ikr_o   = *((real *)((char *)sv + pitch * 40) + threadID_);
        ikr_i   = *((real *)((char *)sv + pitch * 41) + threadID_);
        Jrel_p  = *((real *)((char *)sv + pitch * 42) + threadID_);
    }

    #include "ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.common.c"
}


inline __device__ void RHS_RL_gpu(real *a_, real *b_, real *sv, real *rDY_, real stim_current, real transmurality, real sf_Iks, real *extra_params, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {
    
    // Current modifiers
    real INa_Multiplier = extra_params[0];   
    real INaL_Multiplier = extra_params[1];  
    real INaCa_Multiplier = extra_params[2];  
    real INaK_Multiplier = extra_params[3];  
    real INab_Multiplier = extra_params[4];   
    real Ito_Multiplier = extra_params[5];  
    real IKr_Multiplier = extra_params[6];   
    real IKs_Multiplier = extra_params[7];   
    real IK1_Multiplier = extra_params[8];  
    real IKb_Multiplier = extra_params[9];  
    real IKCa_Multiplier = extra_params[10];  
    real ICaL_Multiplier = extra_params[11];   
    real ICab_Multiplier = extra_params[12];   
    real IpCa_Multiplier = extra_params[13]; 
    real ICaCl_Multiplier = extra_params[14];  
    real IClb_Multiplier = extra_params[15]; 
    real Jrel_Multiplier = extra_params[16]; 
    real Jup_Multiplier = extra_params[17]; 
    real aCaMK_Multiplier = extra_params[18]; 
    real taurelp_Multiplier = extra_params[19];

    // Get the celltype for the current cell
    real celltype = transmurality;
    
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
    real nca;
    real nca_i;
    real ikr_c0;
    real ikr_c1;
    real ikr_c2;
    real ikr_i;
    real ikr_o;
    real xs1;
    real xs2;
    real Jrel_np;
    real Jrel_p;

    if (use_adpt_dt) {
        v       = sv[0];
        nai     = sv[1];
        nass    = sv[2];
        ki      = sv[3];
        kss     = sv[4];
        cai     = sv[5];
        cass    = sv[6];
        cansr   = sv[7];
        cajsr   = sv[8];
        m       = sv[9];
        hp      = sv[10];
        h       = sv[11];
        j       = sv[12];
        jp      = sv[13];
        mL      = sv[14];
        hL      = sv[15];
        hLp     = sv[16];
        a       = sv[17];
        iF      = sv[18];
        iS      = sv[19];
        ap      = sv[20];
        iFp     = sv[21];
        iSp     = sv[22];
        d       = sv[23];
        ff      = sv[24];
        fs      = sv[25];
        fcaf    = sv[26];
        fcas    = sv[27];
        jca     = sv[28];
        nca     = sv[29];
        nca_i   = sv[30];
        ffp     = sv[31];
        fcafp   = sv[32];
        xs1     = sv[33];
        xs2     = sv[34];
        Jrel_np = sv[35];
        CaMKt   = sv[36];
        ikr_c0  = sv[37];
        ikr_c1  = sv[38];
        ikr_c2  = sv[39];
        ikr_o   = sv[40];
        ikr_i   = sv[41];
        Jrel_p  = sv[42];
    } else {
        v       = *((real *)((char *)sv + pitch * 0) + threadID_);
        nai     = *((real *)((char *)sv + pitch * 1) + threadID_);
        nass    = *((real *)((char *)sv + pitch * 2) + threadID_);
        ki      = *((real *)((char *)sv + pitch * 3) + threadID_);
        kss     = *((real *)((char *)sv + pitch * 4) + threadID_);
        cai     = *((real *)((char *)sv + pitch * 5) + threadID_);
        cass    = *((real *)((char *)sv + pitch * 6) + threadID_);
        cansr   = *((real *)((char *)sv + pitch * 7) + threadID_);
        cajsr   = *((real *)((char *)sv + pitch * 8) + threadID_);
        m       = *((real *)((char *)sv + pitch * 9) + threadID_);
        hp      = *((real *)((char *)sv + pitch * 10) + threadID_);
        h       = *((real *)((char *)sv + pitch * 11) + threadID_);
        j       = *((real *)((char *)sv + pitch * 12) + threadID_);
        jp      = *((real *)((char *)sv + pitch * 13) + threadID_);
        mL      = *((real *)((char *)sv + pitch * 14) + threadID_);
        hL      = *((real *)((char *)sv + pitch * 15) + threadID_);
        hLp     = *((real *)((char *)sv + pitch * 16) + threadID_);
        a       = *((real *)((char *)sv + pitch * 17) + threadID_);
        iF      = *((real *)((char *)sv + pitch * 18) + threadID_);
        iS      = *((real *)((char *)sv + pitch * 19) + threadID_);
        ap      = *((real *)((char *)sv + pitch * 20) + threadID_);
        iFp     = *((real *)((char *)sv + pitch * 21) + threadID_);
        iSp     = *((real *)((char *)sv + pitch * 22) + threadID_);
        d       = *((real *)((char *)sv + pitch * 23) + threadID_);
        ff      = *((real *)((char *)sv + pitch * 24) + threadID_);
        fs      = *((real *)((char *)sv + pitch * 25) + threadID_);
        fcaf    = *((real *)((char *)sv + pitch * 26) + threadID_);
        fcas    = *((real *)((char *)sv + pitch * 27) + threadID_);
        jca     = *((real *)((char *)sv + pitch * 28) + threadID_);
        nca     = *((real *)((char *)sv + pitch * 29) + threadID_);
        nca_i   = *((real *)((char *)sv + pitch * 30) + threadID_);
        ffp     = *((real *)((char *)sv + pitch * 31) + threadID_);
        fcafp   = *((real *)((char *)sv + pitch * 32) + threadID_);
        xs1     = *((real *)((char *)sv + pitch * 33) + threadID_);
        xs2     = *((real *)((char *)sv + pitch * 34) + threadID_);
        Jrel_np = *((real *)((char *)sv + pitch * 35) + threadID_);
        CaMKt   = *((real *)((char *)sv + pitch * 36) + threadID_);
        ikr_c0  = *((real *)((char *)sv + pitch * 37) + threadID_);
        ikr_c1  = *((real *)((char *)sv + pitch * 38) + threadID_);
        ikr_c2  = *((real *)((char *)sv + pitch * 39) + threadID_);
        ikr_o   = *((real *)((char *)sv + pitch * 40) + threadID_);
        ikr_i   = *((real *)((char *)sv + pitch * 41) + threadID_);
        Jrel_p  = *((real *)((char *)sv + pitch * 42) + threadID_);
    }

    #include "ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments_RL.common.c"
}
