#include "ToRORd_dynCl_mixed_endo_mid_epi.h"
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
    log_info("Using ToRORd_dynCl_mixed_endo_mid_epi GPU model\n");

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

    // Execution configuration
    const int GRID = (num_volumes + BLOCK_SIZE - 1) / BLOCK_SIZE;

    size_t size = num_volumes * sizeof(real);

    if(use_adpt_h)
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ + 3));
    else
        check_cuda_error(cudaMallocPitch((void **)&(solver->sv), &pitch_h, size, (size_t)NEQ));

    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

    // Get initial condition from extra_data
    real *initial_conditions_endo = NULL;
    real *initial_conditions_epi = NULL;
    real *initial_conditions_mid = NULL;
    real *mapping = NULL;
    real *initial_conditions_endo_device = NULL;
    real *initial_conditions_epi_device = NULL;
    real *initial_conditions_mid_device = NULL;
    real *mapping_device = NULL;

    if(solver->ode_extra_data) 
    {
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
    else
    {
        log_error_and_exit("You must supply a mask function to tag the cells when using this mixed model!\n");
    }

    kernel_set_model_initial_conditions<<<GRID, BLOCK_SIZE>>>(solver->sv,\
                                                            initial_conditions_endo_device, initial_conditions_epi_device, initial_conditions_mid_device,\
                                                            mapping_device, num_volumes);

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

    // Get the mapping array
    uint32_t num_volumes = ode_solver->original_num_cells;
    real *mapping = NULL;
    real *mapping_device = NULL;
    if(ode_solver->ode_extra_data) 
    {
        mapping = (real *)ode_solver->ode_extra_data+NEQ+NEQ+NEQ;
        check_cuda_error(cudaMalloc((void **)&mapping_device, sizeof(real)*num_volumes));
        check_cuda_error(cudaMemcpy(mapping_device, mapping, sizeof(real)*num_volumes, cudaMemcpyHostToDevice));
    }
    else 
    {
        log_error_and_exit("You must supply a mask function to tag the cells when using this mixed model!\n");
    }

    solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, mapping_device);

    check_cuda_error(cudaPeekAtLastError());

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if (mapping_device) check_cuda_error(cudaFree(mapping_device));

}

__global__ void kernel_set_model_initial_conditions(real *sv,\
                                                real *initial_endo, real *initial_epi, real *initial_mid,\
                                                real *mapping, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) 
    {
        for (int i = 0; i < NEQ; i++)
        {
            if (mapping[threadID] == ENDO)
                *((real * )((char *) sv + pitch * i) + threadID) = initial_endo[i];
            else if (mapping[threadID] == EPI)
                *((real * )((char *) sv + pitch * i) + threadID) = initial_epi[i];
            else
                *((real * )((char *) sv + pitch * i) + threadID) = initial_mid[i];
        }
            
        if(use_adpt) 
        {
            *((real *)((char *)sv + pitch * (NEQ)) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * (NEQ+1)) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * (NEQ+2)) + threadID) = 0.0;    // previous dt
        }
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                            uint32_t *cells_to_solve, uint32_t num_cells_to_solve, int num_steps, real *mapping)
{
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

                RHS_gpu(sv, rDY, stim_currents[threadID], mapping[threadID], sv_id, dt);

                // Full Explicit Euler
                for(int i = 0; i < NEQ; i++) {
                    *((real *)((char *)sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
                }

                // Explicit Euler + RushLarsen
                //for(int i = 0; i < 12; i++)
                //    sv[i] = *((real *)((char *)sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
                //for(int i = 12; i < 34; i++)
                //    sv[i] = rDY[i];
                //for(int i = 34; i < 41; i++)
                //    sv[i] = *((real *)((char *)sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
                //for(int i = 41; i < NEQ; i++)
                //    sv[i] = rDY[i];
            }
        } else {
            solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], mapping[threadID], cur_time + max_dt, sv_id);
        }
    }
}

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real mapping, real final_time, int thread_id)
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

    const real __tiny_ = pow(abstol, 2.0f);

    // dt = ((time_new + dt) > final_time) ? (final_time - time_new) : dt;
    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    //#pragma unroll
    for(int i = 0; i < NEQ; i++) {
        sv_local[i] = *((real *)((char *)sv + pitch * i) + thread_id);
    }

    RHS_gpu(sv_local, rDY, stim_curr, mapping, thread_id, dt);
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

        RHS_gpu(sv_local, rDY, stim_curr, mapping, thread_id, dt);
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

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, real mapping, int threadID_, real dt) 
{
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
    real cli;
    real clss;
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
    if (use_adpt)
    {
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
        cli = sv[10];
        clss = sv[11];
        m = sv[12];
        h = sv[13];
        j = sv[14];
        hp = sv[15];
        jp = sv[16];
        mL = sv[17];
        hL = sv[18];
        hLp = sv[19];
        a = sv[20];
        iF = sv[21];
        iS = sv[22];
        ap = sv[23];
        iFp = sv[24];
        iSp = sv[25];
        d = sv[26];
        ff = sv[27];
        fs = sv[28];
        fcaf = sv[29];
        fcas = sv[30];
        jca = sv[31];
        ffp = sv[32];
        fcafp = sv[33];
        nca_ss = sv[34];
        nca_i = sv[35];
        C1 = sv[36];
        C2 = sv[37];
        C3 = sv[38];
        I = sv[39];
        O = sv[40];
        xs1 = sv[41];
        xs2 = sv[42];
        Jrel_np = sv[43];
        Jrel_p = sv[44];
    }
    else
    {
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
        cli = *((real *)((char *)sv + pitch * 10) + threadID_);
        clss = *((real *)((char *)sv + pitch * 11) + threadID_);
        m = *((real *)((char *)sv + pitch * 12) + threadID_);
        h = *((real *)((char *)sv + pitch * 13) + threadID_);
        j = *((real *)((char *)sv + pitch * 14) + threadID_);
        hp = *((real *)((char *)sv + pitch * 15) + threadID_);
        jp = *((real *)((char *)sv + pitch * 16) + threadID_);
        mL = *((real *)((char *)sv + pitch * 17) + threadID_);
        hL = *((real *)((char *)sv + pitch * 18) + threadID_);
        hLp = *((real *)((char *)sv + pitch * 19) + threadID_);
        a = *((real *)((char *)sv + pitch * 20) + threadID_);
        iF = *((real *)((char *)sv + pitch * 21) + threadID_);
        iS = *((real *)((char *)sv + pitch * 22) + threadID_);
        ap = *((real *)((char *)sv + pitch * 23) + threadID_);
        iFp = *((real *)((char *)sv + pitch * 24) + threadID_);
        iSp = *((real *)((char *)sv + pitch * 25) + threadID_);
        d = *((real *)((char *)sv + pitch * 26) + threadID_);
        ff = *((real *)((char *)sv + pitch * 27) + threadID_);
        fs = *((real *)((char *)sv + pitch * 28) + threadID_);
        fcaf = *((real *)((char *)sv + pitch * 29) + threadID_);
        fcas = *((real *)((char *)sv + pitch * 30) + threadID_);
        jca = *((real *)((char *)sv + pitch * 31) + threadID_);
        ffp = *((real *)((char *)sv + pitch * 32) + threadID_);
        fcafp = *((real *)((char *)sv + pitch * 33) + threadID_);
        nca_ss = *((real *)((char *)sv + pitch * 34) + threadID_);
        nca_i = *((real *)((char *)sv + pitch * 35) + threadID_);
        C1 = *((real *)((char *)sv + pitch * 36) + threadID_);
        C2 = *((real *)((char *)sv + pitch * 37) + threadID_);
        C3 = *((real *)((char *)sv + pitch * 38) + threadID_);
        I = *((real *)((char *)sv + pitch * 39) + threadID_);
        O = *((real *)((char *)sv + pitch * 40) + threadID_);
        xs1 = *((real *)((char *)sv + pitch * 41) + threadID_);
        xs2 = *((real *)((char *)sv + pitch * 42) + threadID_);
        Jrel_np = *((real *)((char *)sv + pitch * 43) + threadID_);
        Jrel_p = *((real *)((char *)sv + pitch * 44) + threadID_);
    }

    #include "ToRORd_dynCl_mixed_endo_mid_epi.common.c"
}
