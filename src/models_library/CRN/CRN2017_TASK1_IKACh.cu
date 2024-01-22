#include "../../gpu_utils/gpu_utils.h"
#include <stddef.h>
#include <stdint.h>

#include "CRN2017_TASK1_IKACh.h"

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
    log_info("Using CRN2017 Regions GPU model\n");

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

    real *extra_parameters = NULL;
    real *extra_parameters_device = NULL;

    if(solver->ode_extra_data) {
        extra_parameters = (real *)solver->ode_extra_data;
        check_cuda_error(cudaMalloc((void **)&extra_parameters_device, solver->extra_data_size));
        check_cuda_error(cudaMemcpy(extra_parameters_device, extra_parameters, solver->extra_data_size, cudaMemcpyHostToDevice));
    }
    
    kernel_set_model_initial_conditions <<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes, extra_parameters_device);

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

    // the array cells to solve is passed when we are using and adapative mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) 
    {
        check_cuda_error(cudaMalloc((void **)&cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    // Get the extra_data array
    real *extra_data = NULL;
    real *extra_data_device = NULL;
    if(ode_solver->ode_extra_data) 
    {
        extra_data = (real*)ode_solver->ode_extra_data;
    }
    else 
    {
        log_error_and_exit("You need to specify a mask function when using a mixed model!\n");
    }
    check_cuda_error(cudaMalloc((void **)&extra_data_device, ode_solver->extra_data_size));
    check_cuda_error(cudaMemcpy(extra_data_device, extra_data, ode_solver->extra_data_size, cudaMemcpyHostToDevice));
    
    solve_gpu<<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, extra_data_device);

    check_cuda_error(cudaPeekAtLastError());
    
    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if(extra_data_device) check_cuda_error(cudaFree(extra_data_device));
}

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, real *extra_parameters) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {
    
        if(extra_parameters == NULL)
        {
            *((real * )((char *) sv + pitch * 0) + threadID)  = -82.7618;
            *((real * )((char *) sv + pitch * 1) + threadID)  =   0.0022;
            *((real * )((char *) sv + pitch * 2) + threadID)  =   0.9752;
            *((real * )((char *) sv + pitch * 3) + threadID)  =   0.9844;
            *((real * )((char *) sv + pitch * 4) + threadID)  =   0.0279;
            *((real * )((char *) sv + pitch * 5) + threadID)  =   0.9994;
            *((real * )((char *) sv + pitch * 6) + threadID)  =   1.1219e-4;
            *((real * )((char *) sv + pitch * 7) + threadID)  =   0.9070;
            *((real * )((char *) sv + pitch * 8) + threadID)  =   0.7025;
            *((real * )((char *) sv + pitch * 9) + threadID)  =   0.0042;
            *((real * )((char *) sv + pitch * 10) + threadID) =   0.7816;
            *((real * )((char *) sv + pitch * 11) + threadID) =   0.9488;
            *((real * )((char *) sv + pitch * 12) + threadID) =   0.0044;
            *((real * )((char *) sv + pitch * 13) + threadID) =   0.0184;
            *((real * )((char *) sv + pitch * 14) + threadID) =   0.0000;
            *((real * )((char *) sv + pitch * 15) + threadID) =   0.7719;
            *((real * )((char *) sv + pitch * 16) + threadID) =   1.1355;
            *((real * )((char *) sv + pitch * 17) + threadID) =   0.0000;
            *((real * )((char *) sv + pitch * 18) + threadID) =   1.0000;
            *((real * )((char *) sv + pitch * 19) + threadID) =   0.9993;
            *((real * )((char *) sv + pitch * 20) + threadID) =  11.9137;
            *((real * )((char *) sv + pitch * 21) + threadID) = 140.0000;
            *((real * )((char *) sv + pitch * 22) + threadID) = 138.0694;
            *((real * )((char *) sv + pitch * 23) + threadID) =   5.4000;
            *((real * )((char *) sv + pitch * 24) + threadID) =   1.4810e-4;
            *((real * )((char *) sv + pitch * 25) + threadID) =   1.8000;
        }
        else
        {
            *((real * )((char *) sv + pitch * 0) + threadID)  = extra_parameters[14];
            *((real * )((char *) sv + pitch * 1) + threadID)  = extra_parameters[15];
            *((real * )((char *) sv + pitch * 2) + threadID)  = extra_parameters[16];
            *((real * )((char *) sv + pitch * 3) + threadID)  = extra_parameters[17];
            *((real * )((char *) sv + pitch * 4) + threadID)  = extra_parameters[18];
            *((real * )((char *) sv + pitch * 5) + threadID)  = extra_parameters[19];
            *((real * )((char *) sv + pitch * 6) + threadID)  = extra_parameters[20];
            *((real * )((char *) sv + pitch * 7) + threadID)  = extra_parameters[21];
            *((real * )((char *) sv + pitch * 8) + threadID)  = extra_parameters[22];
            *((real * )((char *) sv + pitch * 9) + threadID)  = extra_parameters[23];
            *((real * )((char *) sv + pitch * 10) + threadID) = extra_parameters[24];
            *((real * )((char *) sv + pitch * 11) + threadID) = extra_parameters[25];
            *((real * )((char *) sv + pitch * 12) + threadID) = extra_parameters[26];
            *((real * )((char *) sv + pitch * 13) + threadID) = extra_parameters[27];
            *((real * )((char *) sv + pitch * 14) + threadID) = extra_parameters[28];
            *((real * )((char *) sv + pitch * 15) + threadID) = extra_parameters[29];
            *((real * )((char *) sv + pitch * 16) + threadID) = extra_parameters[30];
            *((real * )((char *) sv + pitch * 17) + threadID) = extra_parameters[31];
            *((real * )((char *) sv + pitch * 18) + threadID) = extra_parameters[32];
            *((real * )((char *) sv + pitch * 19) + threadID) = extra_parameters[33];
            *((real * )((char *) sv + pitch * 20) + threadID) = extra_parameters[34];
            *((real * )((char *) sv + pitch * 21) + threadID) = extra_parameters[35];
            *((real * )((char *) sv + pitch * 22) + threadID) = extra_parameters[36];
            *((real * )((char *) sv + pitch * 23) + threadID) = extra_parameters[37];
            *((real * )((char *) sv + pitch * 24) + threadID) = extra_parameters[38];
            *((real * )((char *) sv + pitch * 25) + threadID) = extra_parameters[39];
        }

        if(use_adpt) 
        {
            *((real *)((char *)sv + pitch * 26) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * 27) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * 28) + threadID) = 0.0;    // previous dt
        }
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,                         
                            uint32_t *cells_to_solve, uint32_t num_cells_to_solve, int num_steps, real *extra_parameters)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    int offset = 40;
    int mapping = extra_parameters[threadID+offset];

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) 
    {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        if(!use_adpt) 
        {
            real rDY[NEQ];

            for(int n = 0; n < num_steps; ++n) 
            {

                RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, extra_parameters, mapping);

                for(int i = 0; i < NEQ; i++) 
                {
                    *((real *)((char *)sv + pitch * i) + sv_id) =
                        dt * rDY[i] + *((real *)((char *)sv + pitch * i) + sv_id);
                }
            }
        } 
        else 
        {
            solve_forward_euler_gpu_adpt(sv, stim_currents[threadID], cur_time + max_dt, sv_id, extra_parameters, mapping);
        }
    }
}

inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real final_time, int threadID, real *extra_parameters, int mapping) 
{

    #define DT *((real *)((char *)sv + pitch * 26) + threadID)
    #define TIME_NEW *((real *)((char *)sv + pitch * 27) + threadID)
    #define PREVIOUS_DT *((real *)((char *)sv + pitch * 28) + threadID)

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
        sv_local[i] = *((real *)((char *)sv + pitch * i) + threadID);
    }

    RHS_gpu(sv_local, rDY, stim_curr, threadID, dt, extra_parameters, mapping);
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

        RHS_gpu(sv_local, rDY, stim_curr, threadID, dt, extra_parameters, mapping);
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
                    printf("Error: %d: %lf\n", threadID, final_time - time_new);
                    break;
                }
            } else {
                time_new += previous_dt;
            }
        }
    }

    //#pragma unroll
    for(int i = 0; i < NEQ; i++) {
        *((real *)((char *)sv + pitch * i) + threadID) = sv_local[i];
    }

    DT = dt;
    TIME_NEW = time_new;
    PREVIOUS_DT = previous_dt;
}

inline __device__ void RHS_gpu(real *sv, real *rDY, real stim_current, int threadID, real dt, real *extra_parameters, int mapping) 
{
    // State variables
    real Vm;
    real INa_va;
    real INa_vi_1;
    real INa_vi_2;
    real Ito_va;
    real Ito_vi;
    real ICaL_va;
    real ICaL_vi;
    real ICaL_ci;
    real IKur_va;
    real IKur_viS;
    real IKur_viF;
    real IKr_va;
    real IKs_va;
    real IK2P_va;
    real CajSR;
    real CanSR;
    real RyRo;
    real RyRr;
    real RyRi;
    real Nai;
    real Nao;
    real Ki;
    real Ko;
    real Cai;
    real Cao;
    
    if(use_adpt) {
        Vm         = sv[ 0];
        INa_va     = sv[ 1];
        INa_vi_1   = sv[ 2];
        INa_vi_2   = sv[ 3];
        Ito_va     = sv[ 4];
        Ito_vi     = sv[ 5];
        ICaL_va    = sv[ 6];
        ICaL_vi    = sv[ 7];
        ICaL_ci    = sv[ 8];
        IKur_va    = sv[ 9];
        IKur_viS   = sv[10];
        IKur_viF   = sv[11];
        IKr_va     = sv[12];
        IKs_va     = sv[13];
        IK2P_va    = sv[14];
        CajSR      = sv[15];
        CanSR      = sv[16];
        RyRo       = sv[17];
        RyRr       = sv[18];
        RyRi       = sv[19];
        Nai        = sv[20];
        Nao        = sv[21];
        Ki         = sv[22];
        Ko         = sv[23];
        Cai        = sv[24];
        Cao        = sv[25];
    } else {
        Vm         =  *((real*)((char*)sv + pitch * 0) + threadID);
        INa_va     =  *((real*)((char*)sv + pitch * 1) + threadID);
        INa_vi_1   =  *((real*)((char*)sv + pitch * 2) + threadID);
        INa_vi_2   =  *((real*)((char*)sv + pitch * 3) + threadID);
        Ito_va     =  *((real*)((char*)sv + pitch * 4) + threadID);
        Ito_vi     =  *((real*)((char*)sv + pitch * 5) + threadID);
        ICaL_va    =  *((real*)((char*)sv + pitch * 6) + threadID);
        ICaL_vi    =  *((real*)((char*)sv + pitch * 7) + threadID);
        ICaL_ci    =  *((real*)((char*)sv + pitch * 8) + threadID);
        IKur_va    =  *((real*)((char*)sv + pitch * 9) + threadID);
        IKur_viS   =  *((real*)((char*)sv + pitch * 10) + threadID);
        IKur_viF   =  *((real*)((char*)sv + pitch * 11) + threadID);
        IKr_va     =  *((real*)((char*)sv + pitch * 12) + threadID);
        IKs_va     =  *((real*)((char*)sv + pitch * 13) + threadID);
        IK2P_va    =  *((real*)((char*)sv + pitch * 14) + threadID);
        CajSR      =  *((real*)((char*)sv + pitch * 15) + threadID);
        CanSR      =  *((real*)((char*)sv + pitch * 16) + threadID);
        RyRo       =  *((real*)((char*)sv + pitch * 17) + threadID);
        RyRr       =  *((real*)((char*)sv + pitch * 18) + threadID);
        RyRi       =  *((real*)((char*)sv + pitch * 19) + threadID);
        Nai        =  *((real*)((char*)sv + pitch * 20) + threadID);
        Nao        =  *((real*)((char*)sv + pitch * 21) + threadID);
        Ki         =  *((real*)((char*)sv + pitch * 22) + threadID);
        Ko         =  *((real*)((char*)sv + pitch * 23) + threadID);
        Cai        =  *((real*)((char*)sv + pitch * 24) + threadID);
        Cao        =  *((real*)((char*)sv + pitch * 25) + threadID);
    }

    #include "CRN2017_TASK1_IKACh_common.inc.c"
}
