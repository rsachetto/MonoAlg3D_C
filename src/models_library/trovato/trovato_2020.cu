#include "trovato_2020.h"
#include <stddef.h>
#include <stdint.h>

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        // Default initial conditions (1000 beats at BCL=1000ms)
        /*
        *((real * )((char *) sv + pitch * 0) + threadID)      = -86.550102957989600;
        *((real * )((char *) sv + pitch * 1) + threadID)      = 0.005060490773142;
        *((real * )((char *) sv + pitch * 2) + threadID)      = 1.017658486729359e-04;
        *((real * )((char *) sv + pitch * 3) + threadID)      = 8.231857731510893;     
        *((real * )((char *) sv + pitch * 4) + threadID)      = 8.231553454361393;     
        *((real * )((char *) sv + pitch * 5) + threadID)      = 8.231561539013950;
        *((real * )((char *) sv + pitch * 6) + threadID)      = 1.437673447706863e+02;  
        *((real * )((char *) sv + pitch * 7) + threadID)      = 1.437677533510394e+02;
        *((real * )((char *) sv + pitch * 8) + threadID)      = 1.437677771226899e+02;  
        *((real * )((char *) sv + pitch * 9) + threadID)      = 4.360080908582633e-05; 
        *((real * )((char *) sv + pitch * 10) + threadID)      = 1.020101597047671e-04; 
        *((real * )((char *) sv + pitch * 11) + threadID)      = 1.263525645841406;   
        *((real * )((char *) sv + pitch * 12) + threadID)      = 1.248146625349512;    
        *((real * )((char *) sv + pitch * 13) + threadID)      = 1.265185617721750;
        *((real * )((char *) sv + pitch * 14) + threadID)      = 0;                    
        *((real * )((char *) sv + pitch * 15) + threadID)      = 0;
        *((real * )((char *) sv + pitch * 16) + threadID)      = 0.006341207769833;       
        *((real * )((char *) sv + pitch * 17) + threadID)      = 0.788541761218318;       
        *((real * )((char *) sv + pitch * 18) + threadID)      = 0.788474863764949;
        *((real * )((char *) sv + pitch * 19) + threadID)      = 0.790412100577539;       
        *((real * )((char *) sv + pitch * 20) + threadID)      = 0.579594138981772;      
        *((real * )((char *) sv + pitch * 21) + threadID)      = 0.790885181877794;
        *((real * )((char *) sv + pitch * 22) + threadID)      = 0.0;                    
        *((real * )((char *) sv + pitch * 23) + threadID)      = 0.463478975643765;       
        *((real * )((char *) sv + pitch * 24) + threadID)      = 0.240123148015689;
        *((real * )((char *) sv + pitch * 25) + threadID)      = 0.0;                     
        *((real * )((char *) sv + pitch * 26) + threadID)      = 0.649386868151536;        
        *((real * )((char *) sv + pitch * 27) + threadID)      = 0.989963717273401;
        *((real * )((char *) sv + pitch * 28) + threadID)      = 0.0;                     
        *((real * )((char *) sv + pitch * 29) + threadID)      = 0.999999963501871;       
        *((real * )((char *) sv + pitch * 30) + threadID)      = 0.926598485399264;
        *((real * )((char *) sv + pitch * 31) + threadID)      = 0.999999963493016;    
        *((real * )((char *) sv + pitch * 32) + threadID)      = 0.999834100785052;     
        *((real * )((char *) sv + pitch * 33) + threadID)      = 0.999979200703676;
        *((real * )((char *) sv + pitch * 34) + threadID)      = 0.999999963494571;      
        *((real * )((char *) sv + pitch * 35) + threadID)      = 0.999999963612862;
        *((real * )((char *) sv + pitch * 36) + threadID)      = 0.005470852996192;     
        *((real * )((char *) sv + pitch * 37) + threadID)      = 0.0;                     
        *((real * )((char *) sv + pitch * 38) + threadID)      = 0.994211562437775;        
        *((real * )((char *) sv + pitch * 39) + threadID)      = 0.0;
        *((real * )((char *) sv + pitch * 40) + threadID)      = 0.568856244015729;     
        *((real * )((char *) sv + pitch * 41) + threadID)      = 0.191294664752654;      
        *((real * )((char *) sv + pitch * 42) + threadID)      = 0.0;
        *((real * )((char *) sv + pitch * 43) + threadID)      = 0.233014639857230;        
        *((real * )((char *) sv + pitch * 44) + threadID)      = 0.997085416662044;     
        *((real * )((char *) sv + pitch * 45) + threadID)      = 0.466232550741101;
        */
        // Initial conditions (200 beats at BCL=1000ms)
        *((real * )((char *) sv + pitch * 0) + threadID) = -8.668819e+01;
        *((real * )((char *) sv + pitch * 1) + threadID) = 4.962000e-03;
        *((real * )((char *) sv + pitch * 2) + threadID) = 1.020000e-04;
        *((real * )((char *) sv + pitch * 3) + threadID) = 8.240764e+00;
        *((real * )((char *) sv + pitch * 4) + threadID) = 8.240456e+00;
        *((real * )((char *) sv + pitch * 5) + threadID) = 8.240464e+00;
        *((real * )((char *) sv + pitch * 6) + threadID) = 1.437586e+02;
        *((real * )((char *) sv + pitch * 7) + threadID) = 1.437590e+02;
        *((real * )((char *) sv + pitch * 8) + threadID) = 1.437590e+02;
        *((real * )((char *) sv + pitch * 9) + threadID) = 4.400000e-05;
        *((real * )((char *) sv + pitch * 10) + threadID) = 1.020000e-04;
        *((real * )((char *) sv + pitch * 11) + threadID) = 1.260273e+00;
        *((real * )((char *) sv + pitch * 12) + threadID) = 1.245054e+00;
        *((real * )((char *) sv + pitch * 13) + threadID) = 1.261940e+00;
        *((real * )((char *) sv + pitch * 14) + threadID) = 1.530000e-04;
        *((real * )((char *) sv + pitch * 15) + threadID) = 0.000000e+00;
        *((real * )((char *) sv + pitch * 16) + threadID) = 6.321000e-03;
        *((real * )((char *) sv + pitch * 17) + threadID) = 7.887830e-01;
        *((real * )((char *) sv + pitch * 18) + threadID) = 7.887200e-01;
        *((real * )((char *) sv + pitch * 19) + threadID) = 7.906340e-01;
        *((real * )((char *) sv + pitch * 20) + threadID) = 5.799350e-01;
        *((real * )((char *) sv + pitch * 21) + threadID) = 7.911100e-01;
        *((real * )((char *) sv + pitch * 22) + threadID) = 2.420000e-04;
        *((real * )((char *) sv + pitch * 23) + threadID) = 4.638110e-01;
        *((real * )((char *) sv + pitch * 24) + threadID) = 2.404730e-01;
        *((real * )((char *) sv + pitch * 25) + threadID) = 2.730000e-04;
        *((real * )((char *) sv + pitch * 26) + threadID) = 6.497620e-01;
        *((real * )((char *) sv + pitch * 27) + threadID) = 9.899710e-01;
        *((real * )((char *) sv + pitch * 28) + threadID) = 0.000000e+00;
        *((real * )((char *) sv + pitch * 29) + threadID) = 1.000000e+00;
        *((real * )((char *) sv + pitch * 30) + threadID) = 9.268520e-01;
        *((real * )((char *) sv + pitch * 31) + threadID) = 1.000000e+00;
        *((real * )((char *) sv + pitch * 32) + threadID) = 1.000000e+00;
        *((real * )((char *) sv + pitch * 33) + threadID) = 9.999790e-01;
        *((real * )((char *) sv + pitch * 34) + threadID) = 1.000000e+00;
        *((real * )((char *) sv + pitch * 35) + threadID) = 1.000000e+00;
        *((real * )((char *) sv + pitch * 36) + threadID) = 5.431000e-03;
        *((real * )((char *) sv + pitch * 37) + threadID) = 3.040000e-04;
        *((real * )((char *) sv + pitch * 38) + threadID) = 9.942210e-01;
        *((real * )((char *) sv + pitch * 39) + threadID) = 3.290000e-04;
        *((real * )((char *) sv + pitch * 40) + threadID) = 5.707310e-01;
        *((real * )((char *) sv + pitch * 41) + threadID) = 1.909580e-01;
        *((real * )((char *) sv + pitch * 42) + threadID) = 2.230000e-04;
        *((real * )((char *) sv + pitch * 43) + threadID) = 2.334660e-01;
        *((real * )((char *) sv + pitch * 44) + threadID) = 9.970830e-01;
        *((real * )((char *) sv + pitch * 45) + threadID) = 4.655710e-01;

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

    real *stims_currents_device = NULL;
    check_cuda_error(cudaMalloc((void **)&stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    // the array cells to solve is passed when we are using and adaptive mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **)&cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(
            cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }
    
    // Get the extra data array if exists
    uint32_t num_volumes = ode_solver->original_num_cells;
    int num_extra_parameters = 29;
    real extra_par[num_extra_parameters];
    real *extra_par_device = NULL;
    if(ode_solver->ode_extra_data) {
        struct extra_data_for_trovato *extra_data = (struct extra_data_for_trovato*)ode_solver->ode_extra_data;
        extra_par[0]     = extra_data->GNa_Multiplier;
        extra_par[1]     = extra_data->GNaL_Multiplier;
        extra_par[2]     = extra_data->GCaT_Multiplier;
        extra_par[3]     = extra_data->Gto_Multiplier;
        extra_par[4]     = extra_data->Gsus_Multiplier;
        extra_par[5]     = extra_data->Gkr_Multiplier;
        extra_par[6]     = extra_data->Gks_Multiplier;
        extra_par[7]     = extra_data->GfNa_Multiplier;
        extra_par[8]     = extra_data->GfK_Multiplier;
        extra_par[9]     = extra_data->GK1_Multiplier;
        extra_par[10]    = extra_data->GNCX_Multiplier;
        extra_par[11]    = extra_data->GNaK_Multiplier;
        extra_par[12]    = extra_data->INa_Multiplier;
        extra_par[13]    = extra_data->ICaL_Multiplier;
        extra_par[14]    = extra_data->ICaNa_Multiplier;
        extra_par[15]    = extra_data->ICaK_Multiplier;
        extra_par[16]    = extra_data->Ito_Multiplier;
        extra_par[17]    = extra_data->INaL_Multiplier;
        extra_par[18]    = extra_data->IKr_Multiplier;
        extra_par[19]    = extra_data->IKs_Multiplier;
        extra_par[20]    = extra_data->IK1_Multiplier;
        extra_par[21]    = extra_data->INaCa_Multiplier;
        extra_par[22]    = extra_data->INaK_Multiplier;
        extra_par[23]    = extra_data->INab_Multiplier;
        extra_par[24]    = extra_data->ICab_Multiplier;
        extra_par[25]    = extra_data->ICaT_Multiplier;
        extra_par[26]    = extra_data->Isus_Multiplier;
        extra_par[27]    = extra_data->If_Multiplier;
        extra_par[28]    = extra_data->IpCa_Multiplier;
        
        check_cuda_error(cudaMalloc((void **)&extra_par_device, sizeof(real)*num_extra_parameters));
        check_cuda_error(cudaMemcpy(extra_par_device, extra_par, sizeof(real)*num_extra_parameters, cudaMemcpyHostToDevice));
    }
    else {
        for (int i = 0; i < num_extra_parameters; i++) {
            extra_par[i] = 1.0;
        }

        check_cuda_error(cudaMalloc((void **)&extra_par_device, sizeof(real)*num_extra_parameters));
        check_cuda_error(cudaMemcpy(extra_par_device, extra_par, sizeof(real)*num_extra_parameters, cudaMemcpyHostToDevice));
    }

    solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, extra_par_device, num_cells_to_solve,
                                    num_steps, ode_solver->pitch, ode_solver->adaptive, ode_solver->abs_tol,
                                    ode_solver->rel_tol, ode_solver->max_dt);

    check_cuda_error(cudaPeekAtLastError());

    if (stims_currents_device) check_cuda_error(cudaFree(stims_currents_device));
    if (cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if (extra_par_device) check_cuda_error(cudaFree(extra_par_device));
}

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real *extra_params, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt) {

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

    RHS_gpu(sv_local, rDY, stim_curr, extra_params, thread_id, dt, pitch, true);
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

		RHS_gpu(sv_local, rDY, stim_curr, extra_params, thread_id, dt, pitch, true);
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

inline __device__ void solve_rush_larsen_gpu_adpt(real *sv, real stim_curr, real *extra_params, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt) {

    #define DT *((real *)((char *)sv + pitch * (NEQ)) + thread_id)
    #define TIME_NEW *((real *)((char *)sv + pitch * (NEQ+1)) + thread_id)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * (NEQ+2)) + thread_id)

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

    const real __tiny_ = pow(abstol, 2.0);

    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    for(int i = 0; i < NEQ; i++) {
        sv_local[i] = *((real *)((char *)sv + pitch * i) + thread_id);
    }

    RHS_RL_gpu(a_, b_, sv_local, rDY, stim_curr, extra_params, thread_id, dt, pitch, true);
    time_new += dt;

    for(int i = 0; i < NEQ; i++) {
        _k1__[i] = rDY[i];
    }

	while(1) {

		SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(0);        // v        
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(1);        // CaMKt    
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(2);        // cass 
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(3);        // nai  
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(4);        // nasl 
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(5);        // nass 
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(6);        // ki 
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(7);        // kss
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(8);        // ksl
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(9);        // cai
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(10);       // casl
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(11);       // cansr
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(12);       // cajsr
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(13);       // cacsr
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(14);       // Jrel1
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(15);       // Jrel2
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(16);          // m
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(17);          // hf
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(18);          // hs
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(19);          // j
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(20);          // hsp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(21);          // jp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(22);          // mL
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(23);          // hL
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(24);          // hLp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(25);          // a
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(26);          // i1
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(27);          // i2
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(28);          // d
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(29);          // ff
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(30);          // fs
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(31);          // fcaf
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(32);          // fcas
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(33);          // jca
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(34);          // ffp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(35);          // fcafp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(36);       // nca
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(37);          // b
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(38);          // g
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(39);          // xrf
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(40);          // xrs
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(41);          // xs1
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(42);          // xs2
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(43);          // y
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_GPU(44);          // xk1
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_GPU(45);       // u

		time_new += dt;

		RHS_RL_gpu(a_new, b_new, sv_local, rDY, stim_curr, extra_params, thread_id, dt, pitch, true);
		time_new -= dt; // step back

		real greatestError = 0.0, auxError = 0.0;
		real as, bs, f, y_2nd_order;
		SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(0);        // v        
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(1);        // CaMKt    
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(2);        // cass 
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(3);        // nai  
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(4);        // nasl 
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(5);        // nass 
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(6);        // ki 
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(7);        // kss
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(8);        // ksl
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(9);        // cai
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(10);       // casl
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(11);       // cansr
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(12);       // cajsr
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(13);       // cacsr
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(14);       // Jrel1
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(15);       // Jrel2
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(16);          // m
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(17);          // hf
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(18);          // hs
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(19);          // j
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(20);          // hsp
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(21);          // jp
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(22);          // mL
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(23);          // hL
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(24);          // hLp
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(25);          // a
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(26);          // i1
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(27);          // i2
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(28);          // d
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(29);          // ff
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(30);          // fs
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(31);          // fcaf
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(32);          // fcas
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(33);          // jca
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(34);          // ffp
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(35);          // fcafp
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(36);       // nca
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(37);          // b
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(38);          // g
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(39);          // xrf
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(40);          // xrs
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(41);          // xs1
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(42);          // xs2
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(43);          // y
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_GPU(44);          // xk1
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_GPU(45);       // u

		/// adapt the time step
		greatestError += __tiny_;
		previous_dt = dt;

		/// adapt the time step
		dt = dt * sqrt(0.5 * reltol / greatestError);

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

__global__ void solve_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, real *extra_params,
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt,
                          real abstol, real reltol, real max_dt) {
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

                RHS_RL_gpu(a, b, sv, rDY, stim_currents[threadID], extra_params, sv_id, dt, pitch, false);

                // Solve variables based on its type:
                //  Non-linear = Euler
                //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
                SOLVE_EQUATION_EULER_GPU(0);        // v        
                SOLVE_EQUATION_EULER_GPU(1);        // CaMKt    
                SOLVE_EQUATION_EULER_GPU(2);        // cass 
                SOLVE_EQUATION_EULER_GPU(3);        // nai  
                SOLVE_EQUATION_EULER_GPU(4);        // nasl 
                SOLVE_EQUATION_EULER_GPU(5);        // nass 
                SOLVE_EQUATION_EULER_GPU(6);        // ki 
                SOLVE_EQUATION_EULER_GPU(7);        // kss
                SOLVE_EQUATION_EULER_GPU(8);        // ksl
                SOLVE_EQUATION_EULER_GPU(9);        // cai
                SOLVE_EQUATION_EULER_GPU(10);       // casl
                SOLVE_EQUATION_EULER_GPU(11);       // cansr
                SOLVE_EQUATION_EULER_GPU(12);       // cajsr
                SOLVE_EQUATION_EULER_GPU(13);       // cacsr
                SOLVE_EQUATION_EULER_GPU(14);       // Jrel1
                SOLVE_EQUATION_EULER_GPU(15);       // Jrel2
                SOLVE_EQUATION_RUSH_LARSEN_GPU(16); // m
                SOLVE_EQUATION_RUSH_LARSEN_GPU(17); // hf
                SOLVE_EQUATION_RUSH_LARSEN_GPU(18); // hs
                SOLVE_EQUATION_RUSH_LARSEN_GPU(19); // j
                SOLVE_EQUATION_RUSH_LARSEN_GPU(20); // hsp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(21); // jp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(22); // mL
                SOLVE_EQUATION_RUSH_LARSEN_GPU(23); // hL
                SOLVE_EQUATION_RUSH_LARSEN_GPU(24); // hLp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(25); // a
                SOLVE_EQUATION_RUSH_LARSEN_GPU(26); // i1
                SOLVE_EQUATION_RUSH_LARSEN_GPU(27); // i2
                SOLVE_EQUATION_RUSH_LARSEN_GPU(28); // d
                SOLVE_EQUATION_RUSH_LARSEN_GPU(29); // ff
                SOLVE_EQUATION_RUSH_LARSEN_GPU(30); // fs
                SOLVE_EQUATION_RUSH_LARSEN_GPU(31); // fcaf
                SOLVE_EQUATION_RUSH_LARSEN_GPU(32); // fcas
                SOLVE_EQUATION_RUSH_LARSEN_GPU(33); // jca
                SOLVE_EQUATION_RUSH_LARSEN_GPU(34); // ffp
                SOLVE_EQUATION_RUSH_LARSEN_GPU(35); // fcafp
                SOLVE_EQUATION_EULER_GPU(36);       // nca
                SOLVE_EQUATION_RUSH_LARSEN_GPU(37); // b
                SOLVE_EQUATION_RUSH_LARSEN_GPU(38); // g
                SOLVE_EQUATION_RUSH_LARSEN_GPU(39); // xrf
                SOLVE_EQUATION_RUSH_LARSEN_GPU(40); // xrs
                SOLVE_EQUATION_RUSH_LARSEN_GPU(41); // xs1
                SOLVE_EQUATION_RUSH_LARSEN_GPU(42); // xs2
                SOLVE_EQUATION_RUSH_LARSEN_GPU(43); // y
                SOLVE_EQUATION_RUSH_LARSEN_GPU(44); // xk1
                SOLVE_EQUATION_EULER_GPU(45);       // u
            }
        } else {
            //solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], extra_params, cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
            solve_rush_larsen_gpu_adpt(sv, stim_currents[threadID], extra_params, cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt);
        }
    }
}

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, real *extra_params, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {

    // Current modifiers
    real GNa_Multiplier     =  extra_params[0];
    real GNaL_Multiplier    =  extra_params[1];
    real GCaT_Multiplier    =  extra_params[2];
    real Gto_Multiplier     =  extra_params[3];
    real Gsus_Multiplier    =  extra_params[4];
    real Gkr_Multiplier     =  extra_params[5];
    real Gks_Multiplier     =  extra_params[6];
    real GfNa_Multiplier    =  extra_params[7];
    real GfK_Multiplier     =  extra_params[8];
    real GK1_Multiplier     =  extra_params[9];
    real GNCX_Multiplier    =  extra_params[10];
    real GNaK_Multiplier    =  extra_params[11];
    real INa_Multiplier     =  extra_params[12]; 
    real ICaL_Multiplier    =  extra_params[13];
    real ICaNa_Multiplier   =  extra_params[14];
    real ICaK_Multiplier    =  extra_params[15];
    real Ito_Multiplier     =  extra_params[16];
    real INaL_Multiplier    =  extra_params[17];
    real IKr_Multiplier     =  extra_params[18]; 
    real IKs_Multiplier     =  extra_params[19]; 
    real IK1_Multiplier     =  extra_params[20]; 
    real INaCa_Multiplier   =  extra_params[21];
    real INaK_Multiplier    =  extra_params[22];  
    real INab_Multiplier    =  extra_params[23];  
    real ICab_Multiplier    =  extra_params[24];
    real ICaT_Multiplier    =  extra_params[25];
    real Isus_Multiplier    =  extra_params[26];
    real If_Multiplier      =  extra_params[27];
    real IpCa_Multiplier    =  extra_params[28];

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real v;
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
        v = sv[0];
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
        v = *((real *)((char *)sv + pitch * 0) + threadID_);
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

    #include "trovato_2020_common.inc.c"
}

inline __device__ void RHS_RL_gpu(real *a_, real *b_, real *sv, real *rDY_, real stim_current, real *extra_params, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {

    // Current modifiers
    real GNa_Multiplier     =  extra_params[0];
    real GNaL_Multiplier    =  extra_params[1];
    real GCaT_Multiplier    =  extra_params[2];
    real Gto_Multiplier     =  extra_params[3];
    real Gsus_Multiplier    =  extra_params[4];
    real Gkr_Multiplier     =  extra_params[5];
    real Gks_Multiplier     =  extra_params[6];
    real GfNa_Multiplier    =  extra_params[7];
    real GfK_Multiplier     =  extra_params[8];
    real GK1_Multiplier     =  extra_params[9];
    real GNCX_Multiplier    =  extra_params[10];
    real GNaK_Multiplier    =  extra_params[11];
    real INa_Multiplier     =  extra_params[12]; 
    real ICaL_Multiplier    =  extra_params[13];
    real ICaNa_Multiplier   =  extra_params[14];
    real ICaK_Multiplier    =  extra_params[15];
    real Ito_Multiplier     =  extra_params[16];
    real INaL_Multiplier    =  extra_params[17];
    real IKr_Multiplier     =  extra_params[18]; 
    real IKs_Multiplier     =  extra_params[19]; 
    real IK1_Multiplier     =  extra_params[20]; 
    real INaCa_Multiplier   =  extra_params[21];
    real INaK_Multiplier    =  extra_params[22];  
    real INab_Multiplier    =  extra_params[23];  
    real ICab_Multiplier    =  extra_params[24];
    real ICaT_Multiplier    =  extra_params[25];
    real Isus_Multiplier    =  extra_params[26];
    real If_Multiplier      =  extra_params[27];
    real IpCa_Multiplier    =  extra_params[28];

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real v;
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
        v = sv[0];
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
        v = *((real *)((char *)sv + pitch * 0) + threadID_);
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

    #include "trovato_2020_RL_common.inc.c"
}
