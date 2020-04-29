#include "ToRORd_fkatp_endo.h"
#include <stddef.h>
#include <stdint.h>

#define DT *((real *)((char *)sv + pitch * 43) + thread_id)
#define TIME_NEW *((real *)((char *)sv + pitch * 44) + thread_id)
#define PREVIOUS_DT *((real *)((char *)sv + pitch * 45) + thread_id)

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    uint8_t use_adpt_h = (uint8_t)solver->adaptive;

    check_cuda_error(cudaMemcpyToSymbol(use_adpt, &use_adpt_h, sizeof(uint8_t)));
    log_to_stdout_and_file("Using ToRORd_fkatp_endo GPU model\n");

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

extern "C" SOLVE_MODEL_ODES_GPU(solve_model_odes_gpu) {

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

    if(threadID < num_volumes) {

        *((real *)((char *)sv + pitch * 0) + threadID) = -8.876380e+01f; // v millivolt
        *((real *)((char *)sv + pitch * 1) + threadID) = 1.110000e-02f;  // CaMKt millimolar
        *((real *)((char *)sv + pitch * 2) + threadID) = 1.210250e+01f;  // nai millimolar
        *((real *)((char *)sv + pitch * 3) + threadID) = 1.210290e+01f;  // nass millimolar
        *((real *)((char *)sv + pitch * 4) + threadID) = 1.423002e+02f;  // ki millimolar
        *((real *)((char *)sv + pitch * 5) + threadID) = 1.423002e+02f;  // kss millimolar
        *((real *)((char *)sv + pitch * 6) + threadID) = 8.158300e-05f;  // cai millimolar
        *((real *)((char *)sv + pitch * 7) + threadID) = 7.030500e-05f;  // cass millimolar
        *((real *)((char *)sv + pitch * 8) + threadID) = 1.521100e+00f;  // cansr millimolar
        *((real *)((char *)sv + pitch * 9) + threadID) = 1.521400e+00f;  // cajsr millimolar
        *((real *)((char *)sv + pitch * 10) + threadID) = 8.057200e-04f; // m dimensionless
        *((real *)((char *)sv + pitch * 11) + threadID) = 8.286000e-01f; // h dimensionless
        *((real *)((char *)sv + pitch * 12) + threadID) = 8.284000e-01f; // j dimensionless
        *((real *)((char *)sv + pitch * 13) + threadID) = 6.707000e-01f; // hp dimensionless
        *((real *)((char *)sv + pitch * 14) + threadID) = 8.281000e-01f; // jp dimensionless
        *((real *)((char *)sv + pitch * 15) + threadID) = 1.629000e-04f; // mL dimensionless
        *((real *)((char *)sv + pitch * 16) + threadID) = 5.255000e-01f; // hL dimensionless
        *((real *)((char *)sv + pitch * 17) + threadID) = 2.872000e-01f; // hLp dimensionless
        *((real *)((char *)sv + pitch * 18) + threadID) = 9.509800e-04f; // a dimensionless
        *((real *)((char *)sv + pitch * 19) + threadID) = 9.996000e-01f; // iF dimensionless
        *((real *)((char *)sv + pitch * 20) + threadID) = 5.936000e-01f; // iS dimensionless
        *((real *)((char *)sv + pitch * 21) + threadID) = 4.845400e-04f; // ap dimensionless
        *((real *)((char *)sv + pitch * 22) + threadID) = 9.996000e-01f; // iFp dimensionless
        *((real *)((char *)sv + pitch * 23) + threadID) = 6.538000e-01f; // iSp dimensionless
        *((real *)((char *)sv + pitch * 24) + threadID) = 8.108400e-09f; // d dimensionless
        *((real *)((char *)sv + pitch * 25) + threadID) = 1.000000e+00f; // ff dimensionless
        *((real *)((char *)sv + pitch * 26) + threadID) = 9.390000e-01f; // fs dimensionless
        *((real *)((char *)sv + pitch * 27) + threadID) = 1.000000e+00f; // fcaf dimensionless
        *((real *)((char *)sv + pitch * 28) + threadID) = 9.999000e-01f; // fcas dimensionless
        *((real *)((char *)sv + pitch * 29) + threadID) = 1.000000e+00f; // jca dimensionless
        *((real *)((char *)sv + pitch * 30) + threadID) = 1.000000e+00f; // ffp dimensionless
        *((real *)((char *)sv + pitch * 31) + threadID) = 1.000000e+00f; // fcafp dimensionless
        *((real *)((char *)sv + pitch * 32) + threadID) = 6.646200e-04f; // nca_ss dimensionless
        *((real *)((char *)sv + pitch * 33) + threadID) = 1.200000e-03f; // nca_i dimensionless
        *((real *)((char *)sv + pitch * 34) + threadID) = 9.981000e-01f; // C3 dimensionless
        *((real *)((char *)sv + pitch * 35) + threadID) = 8.510900e-04f; // C2 dimensionless
        *((real *)((char *)sv + pitch * 36) + threadID) = 7.034400e-04f; // C1 dimensionless
        *((real *)((char *)sv + pitch * 37) + threadID) = 3.758500e-04f; // O dimensionless
        *((real *)((char *)sv + pitch * 38) + threadID) = 1.328900e-05f; // I dimensionless
        *((real *)((char *)sv + pitch * 39) + threadID) = 2.480000e-01f; // xs1 dimensionless
        *((real *)((char *)sv + pitch * 40) + threadID) = 1.770700e-04f; // xs2 dimensionless
        *((real *)((char *)sv + pitch * 41) + threadID) = 1.612900e-22f; // Jrel_np millimolar_per_millisecond
        *((real *)((char *)sv + pitch * 42) + threadID) = 1.247500e-20f; // Jrel_p millimolar_per_millisecond

        if(use_adpt) {
            *((real *)((char *)sv + pitch * 43) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 44) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 45) + threadID) = 0.0;    // previous dt
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
            solve_Forward_Euler_gpu_adpt(sv, stim_currents[threadID], cur_time + max_dt, sv_id);
        }
    }
}

inline __device__ void solve_Forward_Euler_gpu_adpt(real *sv, real stim_curr, real final_time, int thread_id) {

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
        dt = _beta_safety_ * dt * sqrtf(1.0f / greatestError);

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

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt) {

    // State variables
    real v_old_;
    real CaMKt_old_;
    real nai_old_;
    real nass_old_;
    real ki_old_;
    real kss_old_;
    real cai_old_;
    real cass_old_;
    real cansr_old_;
    real cajsr_old_;
    real m_old_;
    real h_old_;
    real j_old_;
    real hp_old_;
    real jp_old_;
    real mL_old_;
    real hL_old_;
    real hLp_old_;
    real a_old_;
    real iF_old_;
    real iS_old_;
    real ap_old_;
    real iFp_old_;
    real iSp_old_;
    real d_old_;
    real ff_old_;
    real fs_old_;
    real fcaf_old_;
    real fcas_old_;
    real jca_old_;
    real ffp_old_;
    real fcafp_old_;
    real nca_ss_old_;
    real nca_i_old_;
    real C3_old_;
    real C2_old_;
    real C1_old_;
    real O_old_;
    real I_old_;
    real xs1_old_;
    real xs2_old_;
    real Jrel_np_old_;
    real Jrel_p_old_;

    if(use_adpt) {
        v_old_ = sv[0];
        CaMKt_old_ = sv[1];
        nai_old_ = sv[2];
        nass_old_ = sv[3];
        ki_old_ = sv[4];
        kss_old_ = sv[5];
        cai_old_ = sv[6];
        cass_old_ = sv[7];
        cansr_old_ = sv[8];
        cajsr_old_ = sv[9];
        m_old_ = sv[10];
        h_old_ = sv[11];
        j_old_ = sv[12];
        hp_old_ = sv[13];
        jp_old_ = sv[14];
        mL_old_ = sv[15];
        hL_old_ = sv[16];
        hLp_old_ = sv[17];
        a_old_ = sv[18];
        iF_old_ = sv[19];
        iS_old_ = sv[20];
        ap_old_ = sv[21];
        iFp_old_ = sv[22];
        iSp_old_ = sv[23];
        d_old_ = sv[24];
        ff_old_ = sv[25];
        fs_old_ = sv[26];
        fcaf_old_ = sv[27];
        fcas_old_ = sv[28];
        jca_old_ = sv[29];
        ffp_old_ = sv[30];
        fcafp_old_ = sv[31];
        nca_ss_old_ = sv[32];
        nca_i_old_ = sv[33];
        C3_old_ = sv[34];
        C2_old_ = sv[35];
        C1_old_ = sv[36];
        O_old_ = sv[37];
        I_old_ = sv[38];
        xs1_old_ = sv[39];
        xs2_old_ = sv[40];
        Jrel_np_old_ = sv[41];
        Jrel_p_old_ = sv[42];
    } else {
        //    //State variables
        v_old_ = *((real *)((char *)sv + pitch * 0) + threadID_);
        CaMKt_old_ = *((real *)((char *)sv + pitch * 1) + threadID_);
        nai_old_ = *((real *)((char *)sv + pitch * 2) + threadID_);
        nass_old_ = *((real *)((char *)sv + pitch * 3) + threadID_);
        ki_old_ = *((real *)((char *)sv + pitch * 4) + threadID_);
        kss_old_ = *((real *)((char *)sv + pitch * 5) + threadID_);
        cai_old_ = *((real *)((char *)sv + pitch * 6) + threadID_);
        cass_old_ = *((real *)((char *)sv + pitch * 7) + threadID_);
        cansr_old_ = *((real *)((char *)sv + pitch * 8) + threadID_);
        cajsr_old_ = *((real *)((char *)sv + pitch * 9) + threadID_);
        m_old_ = *((real *)((char *)sv + pitch * 10) + threadID_);
        h_old_ = *((real *)((char *)sv + pitch * 11) + threadID_);
        j_old_ = *((real *)((char *)sv + pitch * 12) + threadID_);
        hp_old_ = *((real *)((char *)sv + pitch * 13) + threadID_);
        jp_old_ = *((real *)((char *)sv + pitch * 14) + threadID_);
        mL_old_ = *((real *)((char *)sv + pitch * 15) + threadID_);
        hL_old_ = *((real *)((char *)sv + pitch * 16) + threadID_);
        hLp_old_ = *((real *)((char *)sv + pitch * 17) + threadID_);
        a_old_ = *((real *)((char *)sv + pitch * 18) + threadID_);
        iF_old_ = *((real *)((char *)sv + pitch * 19) + threadID_);
        iS_old_ = *((real *)((char *)sv + pitch * 20) + threadID_);
        ap_old_ = *((real *)((char *)sv + pitch * 21) + threadID_);
        iFp_old_ = *((real *)((char *)sv + pitch * 22) + threadID_);
        iSp_old_ = *((real *)((char *)sv + pitch * 23) + threadID_);
        d_old_ = *((real *)((char *)sv + pitch * 24) + threadID_);
        ff_old_ = *((real *)((char *)sv + pitch * 25) + threadID_);
        fs_old_ = *((real *)((char *)sv + pitch * 26) + threadID_);
        fcaf_old_ = *((real *)((char *)sv + pitch * 27) + threadID_);
        fcas_old_ = *((real *)((char *)sv + pitch * 28) + threadID_);
        jca_old_ = *((real *)((char *)sv + pitch * 29) + threadID_);
        ffp_old_ = *((real *)((char *)sv + pitch * 30) + threadID_);
        fcafp_old_ = *((real *)((char *)sv + pitch * 31) + threadID_);
        nca_ss_old_ = *((real *)((char *)sv + pitch * 32) + threadID_);
        nca_i_old_ = *((real *)((char *)sv + pitch * 33) + threadID_);
        C3_old_ = *((real *)((char *)sv + pitch * 34) + threadID_);
        C2_old_ = *((real *)((char *)sv + pitch * 35) + threadID_);
        C1_old_ = *((real *)((char *)sv + pitch * 36) + threadID_);
        O_old_ = *((real *)((char *)sv + pitch * 37) + threadID_);
        I_old_ = *((real *)((char *)sv + pitch * 38) + threadID_);
        xs1_old_ = *((real *)((char *)sv + pitch * 39) + threadID_);
        xs2_old_ = *((real *)((char *)sv + pitch * 40) + threadID_);
        Jrel_np_old_ = *((real *)((char *)sv + pitch * 41) + threadID_);
        Jrel_p_old_ = *((real *)((char *)sv + pitch * 42) + threadID_);
    }

    #include "ToROrd_common.inc.c"
}