#include "Tworld.h"
#include <stddef.h>
#include <stdint.h>

__constant__  size_t pitch;
__constant__  real abstol;
__constant__  real reltol;
__constant__  real max_dt;
__constant__  real min_dt;
__constant__  uint8_t use_adpt;
size_t pitch_h;

#define sv(i) *((real *)((char *)sv + pitch * (i)) + thread_id)

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt, real *extra_data) {
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
	
	int layer = (int) 			extra_data[thread_id];
    
	if (thread_id < num_volumes) {

        #include "TWorld_human_V9_SS.common.cu"

        if(use_adpt) {
            sv(NEQ) = min_dt; // dt
            sv(NEQ+1) = 0.0;    // time_new
            sv(NEQ+2) = 0.0;    // previous dt
        }
    }
}

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    size_t pitch_h;

    uint8_t use_adpt_dt = (uint8_t)solver->adaptive;

    log_info("Using GPU model implemented in %s\n", __FILE__);

    uint32_t num_volumes = solver->original_num_cells;

    if(use_adpt_dt) {
        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_info("Using Euler model to solve the ODEs\n");
    }

    // execution configuration
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
        check_cuda_error(
            cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
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
	
    solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve,
                                    num_steps, ode_solver->pitch, ode_solver->adaptive, ode_solver->abs_tol,
                                    ode_solver->rel_tol, ode_solver->max_dt, extra_data_device);

    check_cuda_error(cudaPeekAtLastError());

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if(extra_data_device) check_cuda_error(cudaFree(extra_data_device));
}


inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt, 
														int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity) 
{

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

    RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt, pitch, true, 
					layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);
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

		RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt, pitch, true, 
							layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);
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
                          real abstol, real reltol, real max_dt, real *extra_data) {
    const real TOLERANCE = 1e-8;
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;
	
    int layer = (int) 			extra_data[threadID];
    real transmural = (int) 	extra_data[threadID + (1*num_cells_to_solve)];
    real apicobasal = 			extra_data[threadID + (2*num_cells_to_solve)];
    real ischemia = 			extra_data[threadID + (3*num_cells_to_solve)];

	
	int modelVF = (int)	 				extra_data[7*num_cells_to_solve + 0];
	int NZ_severity = (int)				extra_data[7*num_cells_to_solve + 1];
	int ICZ_severity = (int) 			extra_data[7*num_cells_to_solve + 2];

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

                RHS_RL_gpu(a, b, sv, rDY, stim_currents[threadID], sv_id, dt, pitch, false, 
							layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);

                // Solve variables based on its type:
                //  Non-linear = Euler
                //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
                for (int i = 0; i <= 23; i++)
                    SOLVE_EQUATION_EULER_GPU(i);
                for (int i = 24; i <= 39; i++)
                    *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
                for (int i = 40; i <= 44; i++)
                    SOLVE_EQUATION_EULER_GPU(i);
                for (int i = 45; i <= 50; i++)
                    *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
                for (int i = 51; i <= 52; i++)
                    SOLVE_EQUATION_EULER_GPU(i);
                for (int i = 53; i <= 54; i++)
                    *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
                for (int i = 55; i <= 64; i++)
                    SOLVE_EQUATION_EULER_GPU(i);
                for (int i = 65; i <= 67; i++)
                    *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
                for (int i = 68; i <= 74; i++)
                    SOLVE_EQUATION_EULER_GPU(i);
                for (int i = 75; i <= 77; i++)
                    *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
                SOLVE_EQUATION_EULER_GPU(78);
                for (int i = 79; i <= 92; i++)
                    *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }
        } else {
            solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], cur_time + max_dt, sv_id, pitch, abstol,  reltol,  dt,  max_dt, 
								layer, transmural, apicobasal, ischemia,  modelVF, NZ_severity, ICZ_severity);
        }
    }
}

inline __device__ void RHS_RL_gpu (	real *a_, real *b_, real *sv, real *rDY, real stim_current, int thread_id, real dt, size_t pitch, bool use_adpt, 
									int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity) {
    //State variables
    real_cpu v;
    real_cpu ions_na_junc;
    real_cpu ions_na_sl;
    real_cpu ions_na_i;
    real_cpu ions_k_i;
    real_cpu ions_ca_junc;
    real_cpu ions_ca_sl;
    real_cpu ions_ca_i;
    real_cpu ions_cl_i;
    real_cpu ions_ca_SR;
    real_cpu buffers_NaBj;
    real_cpu buffers_NaBsl;
    real_cpu buffers_TnClow;
    real_cpu buffers_TnCHc;
    real_cpu buffers_TnCHm;
    real_cpu buffers_CaM;
    real_cpu buffers_Myosin_ca;
    real_cpu buffers_Myosin_mg;
    real_cpu buffers_SRB;
    real_cpu buffers_SLLj;
    real_cpu buffers_SLLsl;
    real_cpu buffers_SLHj;
    real_cpu buffers_SLHsl;
    real_cpu buffers_Csqn;
    real_cpu m;
    real_cpu h;
    real_cpu j;
    real_cpu mL;
    real_cpu hL;
    real_cpu hLp;
    real_cpu ito_xtos;
    real_cpu ito_ytos;
    real_cpu ito_xtof;
    real_cpu ito_ytof;
    real_cpu ito_xtos_p;
    real_cpu ito_ytos_p;
    real_cpu ito_xtof_p;
    real_cpu ito_ytof_p;
    real_cpu iks_xs_junc;
    real_cpu iks_xs_sl;
    real_cpu c0;
    real_cpu c1;
    real_cpu c2;
    real_cpu o;
    real_cpu i;
    real_cpu d;
    real_cpu ff;
    real_cpu fs;
    real_cpu fcaf;
    real_cpu fcas;
    real_cpu jca;
    real_cpu nca;
    real_cpu nca_i;
    real_cpu ffp;
    real_cpu fcafp;
    real_cpu ical_pureCDI_junc;
    real_cpu ical_pureCDI_sl;
    real_cpu ryr_R;
    real_cpu ryr_O;
    real_cpu ryr_I;
    real_cpu ryr_CaRI;
    real_cpu ryr_R_p;
    real_cpu ryr_O_p;
    real_cpu ryr_I_p;
    real_cpu ryr_CaRI_p;
    real_cpu jrel_icaldep_act;
    real_cpu jrel_icaldep_f1;
    real_cpu jrel_icaldep_f2;
    real_cpu contraction_XS;
    real_cpu contraction_XW;
    real_cpu contraction_Ca_TRPN;
    real_cpu contraction_TmBlocked;
    real_cpu contraction_ZETAS;
    real_cpu contraction_ZETAW;
    real_cpu camk_trap;
    real_cpu camk_f_ICaL;
    real_cpu camk_f_RyR;
    real_cpu camk_f_PLB;
    real_cpu casig_serca_trap;
    real_cpu hp;
    real_cpu jp;
    real_cpu m_P;
    real_cpu h_P;
    real_cpu j_P;
    real_cpu hp_P;
    real_cpu jp_P;
    real_cpu d_P;
    real_cpu ff_P;
    real_cpu fs_P;
    real_cpu fcaf_P;
    real_cpu fcas_P;
    real_cpu fBPf;
    real_cpu fcaBPf;

    if(use_adpt) {
        v = sv[0];
        ions_na_junc = sv[1];
        ions_na_sl = sv[2];
        ions_na_i = sv[3];
        ions_k_i = sv[4];
        ions_ca_junc = sv[5];
        ions_ca_sl = sv[6];
        ions_ca_i = sv[7];
        ions_cl_i = sv[8];
        ions_ca_SR = sv[9];
        buffers_NaBj = sv[10];
        buffers_NaBsl = sv[11];
        buffers_TnClow = sv[12];
        buffers_TnCHc = sv[13];
        buffers_TnCHm = sv[14];
        buffers_CaM = sv[15];
        buffers_Myosin_ca = sv[16];
        buffers_Myosin_mg = sv[17];
        buffers_SRB = sv[18];
        buffers_SLLj = sv[19];
        buffers_SLLsl = sv[20];
        buffers_SLHj = sv[21];
        buffers_SLHsl = sv[22];
        buffers_Csqn = sv[23];
        m = sv[24];
        h = sv[25];
        j = sv[26];
        mL = sv[27];
        hL = sv[28];
        hLp = sv[29];
        ito_xtos = sv[30];
        ito_ytos = sv[31];
        ito_xtof = sv[32];
        ito_ytof = sv[33];
        ito_xtos_p = sv[34];
        ito_ytos_p = sv[35];
        ito_xtof_p = sv[36];
        ito_ytof_p = sv[37];
        iks_xs_junc = sv[38];
        iks_xs_sl = sv[39];
        c0 = sv[40];
        c1 = sv[41];
        c2 = sv[42];
        o = sv[43];
        i = sv[44];
        d = sv[45];
        ff = sv[46];
        fs = sv[47];
        fcaf = sv[48];
        fcas = sv[49];
        jca = sv[50];
        nca = sv[51];
        nca_i = sv[52];
        ffp = sv[53];
        fcafp = sv[54];
        ical_pureCDI_junc = sv[55];
        ical_pureCDI_sl = sv[56];
        ryr_R = sv[57];
        ryr_O = sv[58];
        ryr_I = sv[59];
        ryr_CaRI = sv[60];
        ryr_R_p = sv[61];
        ryr_O_p = sv[62];
        ryr_I_p = sv[63];
        ryr_CaRI_p = sv[64];
        jrel_icaldep_act = sv[65];
        jrel_icaldep_f1 = sv[66];
        jrel_icaldep_f2 = sv[67];
        contraction_XS = sv[68];
        contraction_XW = sv[69];
        contraction_Ca_TRPN = sv[70];
        contraction_TmBlocked = sv[71];
        contraction_ZETAS = sv[72];
        contraction_ZETAW = sv[73];
        camk_trap = sv[74];
        camk_f_ICaL = sv[75];
        camk_f_RyR = sv[76];
        camk_f_PLB = sv[77];
        casig_serca_trap = sv[78];
        hp = sv[79];
        jp = sv[80];
        m_P = sv[81];
        h_P = sv[82];
        j_P = sv[83];
        hp_P = sv[84];
        jp_P = sv[85];
        d_P = sv[86];
        ff_P = sv[87];
        fs_P = sv[88];
        fcaf_P = sv[89];
        fcas_P = sv[90];
        fBPf = sv[91];
        fcaBPf = sv[92];
    } else {
        v = sv(0);
        ions_na_junc = sv(1);
        ions_na_sl = sv(2);
        ions_na_i = sv(3);
        ions_k_i = sv(4);
        ions_ca_junc = sv(5);
        ions_ca_sl = sv(6);
        ions_ca_i = sv(7);
        ions_cl_i = sv(8);
        ions_ca_SR = sv(9);
        buffers_NaBj = sv(10);
        buffers_NaBsl = sv(11);
        buffers_TnClow = sv(12);
        buffers_TnCHc = sv(13);
        buffers_TnCHm = sv(14);
        buffers_CaM = sv(15);
        buffers_Myosin_ca = sv(16);
        buffers_Myosin_mg = sv(17);
        buffers_SRB = sv(18);
        buffers_SLLj = sv(19);
        buffers_SLLsl = sv(20);
        buffers_SLHj = sv(21);
        buffers_SLHsl = sv(22);
        buffers_Csqn = sv(23);
        m = sv(24);
        h = sv(25);
        j = sv(26);
        mL = sv(27);
        hL = sv(28);
        hLp = sv(29);
        ito_xtos = sv(30);
        ito_ytos = sv(31);
        ito_xtof = sv(32);
        ito_ytof = sv(33);
        ito_xtos_p = sv(34);
        ito_ytos_p = sv(35);
        ito_xtof_p = sv(36);
        ito_ytof_p = sv(37);
        iks_xs_junc = sv(38);
        iks_xs_sl = sv(39);
        c0 = sv(40);
        c1 = sv(41);
        c2 = sv(42);
        o = sv(43);
        i = sv(44);
        d = sv(45);
        ff = sv(46);
        fs = sv(47);
        fcaf = sv(48);
        fcas = sv(49);
        jca = sv(50);
        nca = sv(51);
        nca_i = sv(52);
        ffp = sv(53);
        fcafp = sv(54);
        ical_pureCDI_junc = sv(55);
        ical_pureCDI_sl = sv(56);
        ryr_R = sv(57);
        ryr_O = sv(58);
        ryr_I = sv(59);
        ryr_CaRI = sv(60);
        ryr_R_p = sv(61);
        ryr_O_p = sv(62);
        ryr_I_p = sv(63);
        ryr_CaRI_p = sv(64);
        jrel_icaldep_act = sv(65);
        jrel_icaldep_f1 = sv(66);
        jrel_icaldep_f2 = sv(67);
        contraction_XS = sv(68);
        contraction_XW = sv(69);
        contraction_Ca_TRPN = sv(70);
        contraction_TmBlocked = sv(71);
        contraction_ZETAS = sv(72);
        contraction_ZETAW = sv(73);
        camk_trap = sv(74);
        camk_f_ICaL = sv(75);
        camk_f_RyR = sv(76);
        camk_f_PLB = sv(77);
        casig_serca_trap = sv(78);
        hp = sv(79);
        jp = sv(80);
        m_P = sv(81);
        h_P = sv(82);
        j_P = sv(83);
        hp_P = sv(84);
        jp_P = sv(85);
        d_P = sv(86);
        ff_P = sv(87);
        fs_P = sv(88);
        fcaf_P = sv(89);
        fcas_P = sv(90);
        fBPf = sv(91);
        fcaBPf = sv(92);
    }

    #include "Tworld_RL_common.inc.c"
}

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int thread_id, real dt, size_t pitch, bool use_adpt, 
									int layer, real transmural, real apicobasal, real ischemia, int modelVF, int NZ_severity, int ICZ_severity)  {
    //State variables
    real_cpu v;
    real_cpu ions_na_junc;
    real_cpu ions_na_sl;
    real_cpu ions_na_i;
    real_cpu ions_k_i;
    real_cpu ions_ca_junc;
    real_cpu ions_ca_sl;
    real_cpu ions_ca_i;
    real_cpu ions_cl_i;
    real_cpu ions_ca_SR;
    real_cpu buffers_NaBj;
    real_cpu buffers_NaBsl;
    real_cpu buffers_TnClow;
    real_cpu buffers_TnCHc;
    real_cpu buffers_TnCHm;
    real_cpu buffers_CaM;
    real_cpu buffers_Myosin_ca;
    real_cpu buffers_Myosin_mg;
    real_cpu buffers_SRB;
    real_cpu buffers_SLLj;
    real_cpu buffers_SLLsl;
    real_cpu buffers_SLHj;
    real_cpu buffers_SLHsl;
    real_cpu buffers_Csqn;
    real_cpu m;
    real_cpu h;
    real_cpu j;
    real_cpu mL;
    real_cpu hL;
    real_cpu hLp;
    real_cpu ito_xtos;
    real_cpu ito_ytos;
    real_cpu ito_xtof;
    real_cpu ito_ytof;
    real_cpu ito_xtos_p;
    real_cpu ito_ytos_p;
    real_cpu ito_xtof_p;
    real_cpu ito_ytof_p;
    real_cpu iks_xs_junc;
    real_cpu iks_xs_sl;
    real_cpu c0;
    real_cpu c1;
    real_cpu c2;
    real_cpu o;
    real_cpu i;
    real_cpu d;
    real_cpu ff;
    real_cpu fs;
    real_cpu fcaf;
    real_cpu fcas;
    real_cpu jca;
    real_cpu nca;
    real_cpu nca_i;
    real_cpu ffp;
    real_cpu fcafp;
    real_cpu ical_pureCDI_junc;
    real_cpu ical_pureCDI_sl;
    real_cpu ryr_R;
    real_cpu ryr_O;
    real_cpu ryr_I;
    real_cpu ryr_CaRI;
    real_cpu ryr_R_p;
    real_cpu ryr_O_p;
    real_cpu ryr_I_p;
    real_cpu ryr_CaRI_p;
    real_cpu jrel_icaldep_act;
    real_cpu jrel_icaldep_f1;
    real_cpu jrel_icaldep_f2;
    real_cpu contraction_XS;
    real_cpu contraction_XW;
    real_cpu contraction_Ca_TRPN;
    real_cpu contraction_TmBlocked;
    real_cpu contraction_ZETAS;
    real_cpu contraction_ZETAW;
    real_cpu camk_trap;
    real_cpu camk_f_ICaL;
    real_cpu camk_f_RyR;
    real_cpu camk_f_PLB;
    real_cpu casig_serca_trap;
    real_cpu hp;
    real_cpu jp;
    real_cpu m_P;
    real_cpu h_P;
    real_cpu j_P;
    real_cpu hp_P;
    real_cpu jp_P;
    real_cpu d_P;
    real_cpu ff_P;
    real_cpu fs_P;
    real_cpu fcaf_P;
    real_cpu fcas_P;
    real_cpu fBPf;
    real_cpu fcaBPf;

    if(use_adpt) {
        v = sv[0];
        ions_na_junc = sv[1];
        ions_na_sl = sv[2];
        ions_na_i = sv[3];
        ions_k_i = sv[4];
        ions_ca_junc = sv[5];
        ions_ca_sl = sv[6];
        ions_ca_i = sv[7];
        ions_cl_i = sv[8];
        ions_ca_SR = sv[9];
        buffers_NaBj = sv[10];
        buffers_NaBsl = sv[11];
        buffers_TnClow = sv[12];
        buffers_TnCHc = sv[13];
        buffers_TnCHm = sv[14];
        buffers_CaM = sv[15];
        buffers_Myosin_ca = sv[16];
        buffers_Myosin_mg = sv[17];
        buffers_SRB = sv[18];
        buffers_SLLj = sv[19];
        buffers_SLLsl = sv[20];
        buffers_SLHj = sv[21];
        buffers_SLHsl = sv[22];
        buffers_Csqn = sv[23];
        m = sv[24];
        h = sv[25];
        j = sv[26];
        mL = sv[27];
        hL = sv[28];
        hLp = sv[29];
        ito_xtos = sv[30];
        ito_ytos = sv[31];
        ito_xtof = sv[32];
        ito_ytof = sv[33];
        ito_xtos_p = sv[34];
        ito_ytos_p = sv[35];
        ito_xtof_p = sv[36];
        ito_ytof_p = sv[37];
        iks_xs_junc = sv[38];
        iks_xs_sl = sv[39];
        c0 = sv[40];
        c1 = sv[41];
        c2 = sv[42];
        o = sv[43];
        i = sv[44];
        d = sv[45];
        ff = sv[46];
        fs = sv[47];
        fcaf = sv[48];
        fcas = sv[49];
        jca = sv[50];
        nca = sv[51];
        nca_i = sv[52];
        ffp = sv[53];
        fcafp = sv[54];
        ical_pureCDI_junc = sv[55];
        ical_pureCDI_sl = sv[56];
        ryr_R = sv[57];
        ryr_O = sv[58];
        ryr_I = sv[59];
        ryr_CaRI = sv[60];
        ryr_R_p = sv[61];
        ryr_O_p = sv[62];
        ryr_I_p = sv[63];
        ryr_CaRI_p = sv[64];
        jrel_icaldep_act = sv[65];
        jrel_icaldep_f1 = sv[66];
        jrel_icaldep_f2 = sv[67];
        contraction_XS = sv[68];
        contraction_XW = sv[69];
        contraction_Ca_TRPN = sv[70];
        contraction_TmBlocked = sv[71];
        contraction_ZETAS = sv[72];
        contraction_ZETAW = sv[73];
        camk_trap = sv[74];
        camk_f_ICaL = sv[75];
        camk_f_RyR = sv[76];
        camk_f_PLB = sv[77];
        casig_serca_trap = sv[78];
        hp = sv[79];
        jp = sv[80];
        m_P = sv[81];
        h_P = sv[82];
        j_P = sv[83];
        hp_P = sv[84];
        jp_P = sv[85];
        d_P = sv[86];
        ff_P = sv[87];
        fs_P = sv[88];
        fcaf_P = sv[89];
        fcas_P = sv[90];
        fBPf = sv[91];
        fcaBPf = sv[92];
    } else {
        v = sv(0);
        ions_na_junc = sv(1);
        ions_na_sl = sv(2);
        ions_na_i = sv(3);
        ions_k_i = sv(4);
        ions_ca_junc = sv(5);
        ions_ca_sl = sv(6);
        ions_ca_i = sv(7);
        ions_cl_i = sv(8);
        ions_ca_SR = sv(9);
        buffers_NaBj = sv(10);
        buffers_NaBsl = sv(11);
        buffers_TnClow = sv(12);
        buffers_TnCHc = sv(13);
        buffers_TnCHm = sv(14);
        buffers_CaM = sv(15);
        buffers_Myosin_ca = sv(16);
        buffers_Myosin_mg = sv(17);
        buffers_SRB = sv(18);
        buffers_SLLj = sv(19);
        buffers_SLLsl = sv(20);
        buffers_SLHj = sv(21);
        buffers_SLHsl = sv(22);
        buffers_Csqn = sv(23);
        m = sv(24);
        h = sv(25);
        j = sv(26);
        mL = sv(27);
        hL = sv(28);
        hLp = sv(29);
        ito_xtos = sv(30);
        ito_ytos = sv(31);
        ito_xtof = sv(32);
        ito_ytof = sv(33);
        ito_xtos_p = sv(34);
        ito_ytos_p = sv(35);
        ito_xtof_p = sv(36);
        ito_ytof_p = sv(37);
        iks_xs_junc = sv(38);
        iks_xs_sl = sv(39);
        c0 = sv(40);
        c1 = sv(41);
        c2 = sv(42);
        o = sv(43);
        i = sv(44);
        d = sv(45);
        ff = sv(46);
        fs = sv(47);
        fcaf = sv(48);
        fcas = sv(49);
        jca = sv(50);
        nca = sv(51);
        nca_i = sv(52);
        ffp = sv(53);
        fcafp = sv(54);
        ical_pureCDI_junc = sv(55);
        ical_pureCDI_sl = sv(56);
        ryr_R = sv(57);
        ryr_O = sv(58);
        ryr_I = sv(59);
        ryr_CaRI = sv(60);
        ryr_R_p = sv(61);
        ryr_O_p = sv(62);
        ryr_I_p = sv(63);
        ryr_CaRI_p = sv(64);
        jrel_icaldep_act = sv(65);
        jrel_icaldep_f1 = sv(66);
        jrel_icaldep_f2 = sv(67);
        contraction_XS = sv(68);
        contraction_XW = sv(69);
        contraction_Ca_TRPN = sv(70);
        contraction_TmBlocked = sv(71);
        contraction_ZETAS = sv(72);
        contraction_ZETAW = sv(73);
        camk_trap = sv(74);
        camk_f_ICaL = sv(75);
        camk_f_RyR = sv(76);
        camk_f_PLB = sv(77);
        casig_serca_trap = sv(78);
        hp = sv(79);
        jp = sv(80);
        m_P = sv(81);
        h_P = sv(82);
        j_P = sv(83);
        hp_P = sv(84);
        jp_P = sv(85);
        d_P = sv(86);
        ff_P = sv(87);
        fs_P = sv(88);
        fcaf_P = sv(89);
        fcas_P = sv(90);
        fBPf = sv(91);
        fcaBPf = sv(92);
    }

    #include "Tworld_common.inc.c"
}
