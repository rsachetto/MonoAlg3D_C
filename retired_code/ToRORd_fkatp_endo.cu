#include "ToRORd_fkatp_endo.h"
#include <stddef.h>
#include <stdint.h>
#include "../gpu_utils/gpu_utils.h"

#define DT *((real*)((char*)sv + pitch * 43) + thread_id)
#define TIME_NEW *((real*)((char*)sv + pitch * 44) + thread_id)
#define PREVIOUS_DT  *((real*)((char*)sv + pitch * 45) + thread_id)

__constant__  int matlab;

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    log_to_stdout_and_file("Using ToRORd_fkatp_endo GPU model\n");

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(solver->sv), &pitch_h, size, (size_t )NEQ + 3));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

    int matlab_h = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, matlab_h, ode_extra_config, "matlab");

    check_cuda_error(cudaMemcpyToSymbol(matlab, &matlab_h, sizeof(int)));

    if(matlab_h == 0)
        kernel_set_model_initial_conditions <<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes);
    else
        kernel_set_model_initial_conditions_matlab <<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_gpu) {

    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));


    //the array cells to solve is passed when we are using and adapative mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    solve_gpu <<<GRID, BLOCK_SIZE>>>(current_t, dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_initial_conditions_matlab(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

         *((real * )((char *) sv + pitch * 0) + threadID)  = -88.8691566357934;     //v
         *((real * )((char *) sv + pitch * 1) + threadID)  = 12.0996647655188;      //nai
         *((real * )((char *) sv + pitch * 2) + threadID)  = 12.1000028563765;      //nass
         *((real * )((char *) sv + pitch * 3) + threadID)  = 142.412524737626;      //ki
         *((real * )((char *) sv + pitch * 4) + threadID)  = 142.412481425842;      //kss
         *((real * )((char *) sv + pitch * 5) + threadID)  = 7.45541572746214e-05;  //cai
         *((real * )((char *) sv + pitch * 6) + threadID)  = 6.50418928341426e-05;  //cass
         *((real * )((char *) sv + pitch * 7) + threadID)  = 1.53037019085812;      //cansr
         *((real * )((char *) sv + pitch * 8) + threadID)  = 1.52803094224238;      //cajsr
         *((real * )((char *) sv + pitch * 9) + threadID)  = 0.000787657400526199;  //m
         *((real * )((char *) sv + pitch * 10) + threadID) = 0.674096901201792;     //hp
         *((real * )((char *) sv + pitch * 11) + threadID) = 0.830658198588696;     //h
         *((real * )((char *) sv + pitch * 12) + threadID) = 0.830466744399495;     //j
         *((real * )((char *) sv + pitch * 13) + threadID) = 0.830093612199637;     //jp
         *((real * )((char *) sv + pitch * 14) + threadID) = 0.000159670117055769;  //mL
         *((real * )((char *) sv + pitch * 15) + threadID) = 0.528261721740178;     //hL
         *((real * )((char *) sv + pitch * 16) + threadID) = 0.288775833197764;     //hLp
         *((real * )((char *) sv + pitch * 17) + threadID) = 0.000944249645410894;  //a
         *((real * )((char *) sv + pitch * 18) + threadID) = 0.999616956857814;     //iF
         *((real * )((char *) sv + pitch * 19) + threadID) = 0.593680589620082;     //iS
         *((real * )((char *) sv + pitch * 20) + threadID) = 0.000481107253796778;  //ap
         *((real * )((char *) sv + pitch * 21) + threadID) = 0.999616964658062;     //iFp
         *((real * )((char *) sv + pitch * 22) + threadID) = 0.654092074678260;     //iSp
         *((real * )((char *) sv + pitch * 23) + threadID) = 8.86091322819384e-29;  //d
         *((real * )((char *) sv + pitch * 24) + threadID) = 0.999999992783113;     //ff
         *((real * )((char *) sv + pitch * 25) + threadID) = 0.938965241412012;     //fs
         *((real * )((char *) sv + pitch * 26) + threadID) = 0.999999992783179;     //fcaf
         *((real * )((char *) sv + pitch * 27) + threadID) = 0.999900458262832;     //fcas
         *((real * )((char *) sv + pitch * 28) + threadID) = 0.999977476316330;     //jca
         *((real * )((char *) sv + pitch * 29) + threadID) = 0.000492094765239740;  //nca
         *((real * )((char *) sv + pitch * 30) + threadID) = 0.000833711885764158;  //nca_i
         *((real * )((char *) sv + pitch * 31) + threadID) = 0.999999992566681;     //ffp
         *((real * )((char *) sv + pitch * 32) + threadID) = 0.999999992766279;     //fcafp
         *((real * )((char *) sv + pitch * 33) + threadID) = 0.247156543918935;     //xs1
         *((real * )((char *) sv + pitch * 34) + threadID) = 0.000175017075236424;  //xs2
         *((real * )((char *) sv + pitch * 35) + threadID) = 3.90843796133124e-24;  //Jrel_np
         *((real * )((char *) sv + pitch * 36) + threadID) = 0.0110752904836162;    //CaMKt
         *((real * )((char *) sv + pitch * 37) + threadID) = 0.998073652444028;     //ikr_c0
         *((real * )((char *) sv + pitch * 38) + threadID) = 0.000844745297078649;  //ikr_c1
         *((real * )((char *) sv + pitch * 39) + threadID) = 0.000698171876592920;  //ikr_c2
         *((real * )((char *) sv + pitch * 40) + threadID) = 0.000370404872169913;  //ikr_o
         *((real * )((char *) sv + pitch * 41) + threadID) = 1.30239063420973e-05;  //ikr_i
         *((real * )((char *) sv + pitch * 42) + threadID) = -1.88428892080206e-22; //Jrel_p

        *((real*)((char*)sv + pitch * 43) + threadID) = 0.05;	 // dt
        *((real*)((char*)sv + pitch * 44) + threadID) = 0.0;	 // time_new
        *((real*)((char*)sv + pitch * 45) + threadID) = 0.0;	 // previous dt
    }
}

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

//           *((real * )((char *) sv + pitch * 0) + threadID) = -88.8691566357934; //v millivolt
//           *((real * )((char *) sv + pitch * 1) + threadID) = 0.0110752904836162; //CaMKt millimolar
//           *((real * )((char *) sv + pitch * 2) + threadID) = 12.0996647655188; //nai millimolar
//           *((real * )((char *) sv + pitch * 3) + threadID) = 12.1000028563765; //nass millimolar
//           *((real * )((char *) sv + pitch * 4) + threadID) = 142.412524737626; //ki millimolar
//           *((real * )((char *) sv + pitch * 5) + threadID) = 142.412481425842; //kss millimolar
//           *((real * )((char *) sv + pitch * 6) + threadID) = 7.45541572746214e-05; //cai millimolar
//           *((real * )((char *) sv + pitch * 7) + threadID) = 6.50418928341426e-05; //cass millimolar
//           *((real * )((char *) sv + pitch * 8) + threadID) = 1.53037019085812; //cansr millimolar
//           *((real * )((char *) sv + pitch * 9) + threadID) = 1.52803094224238; //cajsr millimolar
//           *((real * )((char *) sv + pitch * 10) + threadID) = 0.000787657400526199; //m dimensionless
//           *((real * )((char *) sv + pitch * 11) + threadID) = 0.830658198588696; //h dimensionless
//           *((real * )((char *) sv + pitch * 12) + threadID) = 0.830466744399495; //j dimensionless
//           *((real * )((char *) sv + pitch * 13) + threadID) = 0.674096901201792; //hp dimensionless
//           *((real * )((char *) sv + pitch * 14) + threadID) = 0.830093612199637; //jp dimensionless
//           *((real * )((char *) sv + pitch * 15) + threadID) = 0.000159670117055769; //mL dimensionless
//           *((real * )((char *) sv + pitch * 16) + threadID) = 0.528261721740178; //hL dimensionless
//           *((real * )((char *) sv + pitch * 17) + threadID) = 0.288775833197764; //hLp dimensionless
//           *((real * )((char *) sv + pitch * 18) + threadID) = 0.000944249645410894; //a dimensionless
//           *((real * )((char *) sv + pitch * 19) + threadID) = 0.999616956857814;  //iF dimensionless
//           *((real * )((char *) sv + pitch * 20) + threadID) = 0.593680589620082; //iS dimensionless
//           *((real * )((char *) sv + pitch * 21) + threadID) = 0.000481107253796778; //ap dimensionless
//           *((real * )((char *) sv + pitch * 22) + threadID) = 0.999616964658062; //iFp dimensionless
//           *((real * )((char *) sv + pitch * 23) + threadID) = 0.654092074678260; //iSp dimensionless
//           *((real * )((char *) sv + pitch * 24) + threadID) = 8.86091322819384e-29; //d dimensionless //different in model and matlab
//           *((real * )((char *) sv + pitch * 25) + threadID) = 0.999999992783113; //ff dimensionless
//           *((real * )((char *) sv + pitch * 26) + threadID) = 0.938965241412012; //fs dimensionless
//           *((real * )((char *) sv + pitch * 27) + threadID) = 0.999999992783179; //fcaf dimensionless
//           *((real * )((char *) sv + pitch * 28) + threadID) = 0.999900458262832; //fcas dimensionless
//           *((real * )((char *) sv + pitch * 29) + threadID) = 0.999977476316330; //jca dimensionless
//           *((real * )((char *) sv + pitch * 30) + threadID) = 0.999999992566681; //ffp dimensionless
//           *((real * )((char *) sv + pitch * 31) + threadID) = 0.999999992766279; //fcafp dimensionless
//           *((real * )((char *) sv + pitch * 32) + threadID) = 0.000492094765239740; //nca_ss dimensionless //different in model and matlab
//           *((real * )((char *) sv + pitch * 33) + threadID) = 0.000833711885764158;; //nca_i dimensionless //different in model and matlab
//           *((real * )((char *) sv + pitch * 34) + threadID) = 0.998073652444028; //C3 dimensionless
//           *((real * )((char *) sv + pitch * 35) + threadID) = 0.000844745297078649; //C2 dimensionless
//           *((real * )((char *) sv + pitch * 36) + threadID) = 0.000698171876592920; //C1 dimensionless
//           *((real * )((char *) sv + pitch * 37) + threadID) = 0.000370404872169913; //O dimensionless
//           *((real * )((char *) sv + pitch * 38) + threadID) = 1.30239063420973e-05; //I dimensionless
//           *((real * )((char *) sv + pitch * 39) + threadID) = 0.247156543918935; //xs1 dimensionless
//           *((real * )((char *) sv + pitch * 40) + threadID) = 0.000175017075236424; //xs2 dimensionless
//           *((real * )((char *) sv + pitch * 41) + threadID) = 3.90843796133124e-24; //Jrel_np millimolar_per_millisecond //different in model and matlab
//           *((real * )((char *) sv + pitch * 42) + threadID) = -1.88428892080206e-22; //Jrel_p millimolar_per_millisecond //different in model and matlab

        *((real * )((char *) sv + pitch * 0) + threadID) = -8.876380e+01f; //v millivolt
        *((real * )((char *) sv + pitch * 1) + threadID) = 1.110000e-02f; //CaMKt millimolar
        *((real * )((char *) sv + pitch * 2) + threadID) = 1.210250e+01f; //nai millimolar
        *((real * )((char *) sv + pitch * 3) + threadID) = 1.210290e+01f; //nass millimolar
        *((real * )((char *) sv + pitch * 4) + threadID) = 1.423002e+02f; //ki millimolar
        *((real * )((char *) sv + pitch * 5) + threadID) = 1.423002e+02f; //kss millimolar
        *((real * )((char *) sv + pitch * 6) + threadID) = 8.158300e-05f; //cai millimolar
        *((real * )((char *) sv + pitch * 7) + threadID) = 7.030500e-05f; //cass millimolar
        *((real * )((char *) sv + pitch * 8) + threadID) = 1.521100e+00f; //cansr millimolar
        *((real * )((char *) sv + pitch * 9) + threadID) = 1.521400e+00f; //cajsr millimolar
        *((real * )((char *) sv + pitch * 10) + threadID) = 8.057200e-04f; //m dimensionless
        *((real * )((char *) sv + pitch * 11) + threadID) = 8.286000e-01f; //h dimensionless
        *((real * )((char *) sv + pitch * 12) + threadID) = 8.284000e-01f; //j dimensionless
        *((real * )((char *) sv + pitch * 13) + threadID) = 6.707000e-01f; //hp dimensionless
        *((real * )((char *) sv + pitch * 14) + threadID) = 8.281000e-01f; //jp dimensionless
        *((real * )((char *) sv + pitch * 15) + threadID) = 1.629000e-04f; //mL dimensionless
        *((real * )((char *) sv + pitch * 16) + threadID) = 5.255000e-01f; //hL dimensionless
        *((real * )((char *) sv + pitch * 17) + threadID) = 2.872000e-01f; //hLp dimensionless
        *((real * )((char *) sv + pitch * 18) + threadID) = 9.509800e-04f; //a dimensionless
        *((real * )((char *) sv + pitch * 19) + threadID) = 9.996000e-01f; //iF dimensionless
        *((real * )((char *) sv + pitch * 20) + threadID) = 5.936000e-01f; //iS dimensionless
        *((real * )((char *) sv + pitch * 21) + threadID) = 4.845400e-04f; //ap dimensionless
        *((real * )((char *) sv + pitch * 22) + threadID) = 9.996000e-01f; //iFp dimensionless
        *((real * )((char *) sv + pitch * 23) + threadID) = 6.538000e-01f; //iSp dimensionless
        *((real * )((char *) sv + pitch * 24) + threadID) = 8.108400e-09f; //d dimensionless
        *((real * )((char *) sv + pitch * 25) + threadID) = 1.000000e+00f; //ff dimensionless
        *((real * )((char *) sv + pitch * 26) + threadID) = 9.390000e-01f; //fs dimensionless
        *((real * )((char *) sv + pitch * 27) + threadID) = 1.000000e+00f; //fcaf dimensionless
        *((real * )((char *) sv + pitch * 28) + threadID) = 9.999000e-01f; //fcas dimensionless
        *((real * )((char *) sv + pitch * 29) + threadID) = 1.000000e+00f; //jca dimensionless
        *((real * )((char *) sv + pitch * 30) + threadID) = 1.000000e+00f; //ffp dimensionless
        *((real * )((char *) sv + pitch * 31) + threadID) = 1.000000e+00f; //fcafp dimensionless
        *((real * )((char *) sv + pitch * 32) + threadID) = 6.646200e-04f; //nca_ss dimensionless
        *((real * )((char *) sv + pitch * 33) + threadID) = 1.200000e-03f; //nca_i dimensionless
        *((real * )((char *) sv + pitch * 34) + threadID) = 9.981000e-01f; //C3 dimensionless
        *((real * )((char *) sv + pitch * 35) + threadID) = 8.510900e-04f; //C2 dimensionless
        *((real * )((char *) sv + pitch * 36) + threadID) = 7.034400e-04f; //C1 dimensionless
        *((real * )((char *) sv + pitch * 37) + threadID) = 3.758500e-04f; //O dimensionless
        *((real * )((char *) sv + pitch * 38) + threadID) = 1.328900e-05f; //I dimensionless
        *((real * )((char *) sv + pitch * 39) + threadID) = 2.480000e-01f; //xs1 dimensionless
        *((real * )((char *) sv + pitch * 40) + threadID) = 1.770700e-04f; //xs2 dimensionless
        *((real * )((char *) sv + pitch * 41) + threadID) = 1.612900e-22f; //Jrel_np millimolar_per_millisecond
        *((real * )((char *) sv + pitch * 42) + threadID) = 1.247500e-20f; //Jrel_p millimolar_per_millisecond

        *((real*)((char*)sv + pitch * 43) + threadID) = 0.05;	 // dt
        *((real*)((char*)sv + pitch * 44) + threadID) = 0.0;	 // time_new
        *((real*)((char*)sv + pitch * 45) + threadID) = 0.0;	 // previous dt
    }
}

inline __device__ void solve_Forward_Euler_gpu_adpt(real *sv, real stim_curr, real final_time, real min_step, real max_step, int thread_id);
// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real cur_time, real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps)
{
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell model
    if(threadID < num_cells_to_solve) {
        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

//        real rDY[NEQ];

//        for (int n = 0; n < num_steps; ++n) {
//
//            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt);
//
//            for(int i = 0; i < NEQ; i++) {
//                *((real *) ((char *) sv + pitch * i) + sv_id) = dt*rDY[i] + *((real *) ((char *) sv + pitch * i) + sv_id);
//            }
//
//        }

        solve_Forward_Euler_gpu_adpt(sv, stim_currents[threadID], cur_time + dt, 0.0001, 0.05, sv_id);

    }
}

inline __device__ void solve_Forward_Euler_gpu_adpt(real *sv, real stim_curr, real final_time, real min_step, real max_step, int thread_id) {

    real rDY[NEQ];

    real _tolerances_[NEQ];
    real _aux_tol = 0.0;
    real dt = DT;
    real time_new = TIME_NEW;
    real previous_dt  = PREVIOUS_DT;

    real edos_old_aux_[NEQ];
    real edos_new_euler_[NEQ];
    real _k1__[NEQ];
    real _k2__ [NEQ];
    real _k_aux__[NEQ];
    real sv_local[NEQ];

    const real _beta_safety_ = 0.8;
    const real reltol = 1e-5;
    const real abstol = 1e-5;
    const real __tiny_ = powf(abstol, 2.0f);

    //dt = ((time_new + dt) > final_time) ? (final_time - time_new) : dt;
    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

//#pragma unroll
    for(int i = 0; i < NEQ; i++) {
        sv_local[i] = *((real*)((char*)sv + pitch * i) + thread_id);
    }

    if(matlab == 0)
        RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt);
    else
        RHS_gpu_matlab(sv_local, rDY, stim_curr, thread_id, dt);
    time_new += dt;

//#pragma unroll
    for(int i = 0; i < NEQ; i++) {
        _k1__[i] = rDY[i];
    }

    int count = 0;

    int count_limit = (final_time-time_new)/min_step;

    int aux_count_limit =  count_limit + 2000000;

    if(aux_count_limit > 0) {
        count_limit = aux_count_limit;
    }

    while (1) {

        for (int i = 0; i < NEQ; i++) {
            //stores the old variables in a vector
            edos_old_aux_[i] = sv_local[i];
            // //computes euler method
            edos_new_euler_[i] = _k1__[i] * dt + edos_old_aux_[i];
            //steps ahead to compute the rk2 method
            sv_local[i] = edos_new_euler_[i];
        }

        time_new += dt;

        if(matlab == 0)
            RHS_gpu(sv_local, rDY, stim_curr, thread_id, dt);
        else
            RHS_gpu_matlab(sv_local, rDY, stim_curr, thread_id, dt);
        time_new -= dt; //step back

        real greatestError = 0.0, auxError = 0.0;
//#pragma unroll
        for (int i = 0; i < NEQ; i++) {

            //stores the new evaluation
            _k2__[i] = rDY[i];
            _aux_tol = fabs(edos_new_euler_[i]) * reltol;
            _tolerances_[i] = (abstol > _aux_tol) ? abstol : _aux_tol;

            //finds the greatest error between  the steps
            auxError = fabs(((dt / 2.0) * (_k1__[i] - _k2__[i])) / _tolerances_[i]);

            greatestError = (auxError > greatestError) ? auxError : greatestError;

        }

        ///adapt the time step
        greatestError += __tiny_;
        previous_dt = dt;
        ///adapt the time step
        dt = _beta_safety_ * dt * sqrtf(1.0f / greatestError);

        if (time_new + dt > final_time) {
            dt = final_time - time_new;
        }

        //it doesn't accept the solution
        if ( count < count_limit  && (greatestError >= 1.0f)) {
            // restore the old values to do it again
            for (int i = 0; i < NEQ; i++) {
                sv_local[i] = edos_old_aux_[i];
            }
            count++;
            //throw the results away and compute again
        }
        else {
            //count = 0;

            //if(greatestError >=1.0) {
            //    printf("Thread //d,accepting solution with error > //lf \n", threadID, greatestError);
            //}

            //it accepts the solutions
            //int aux = (dt > max_step && max_step != 0);
            //dt = (aux) ? max_step : dt;

            if (dt < min_step) {
                dt = min_step;
            }

            else if (dt > max_step && max_step != 0) {
                dt = max_step;
            }

            if (time_new + dt > final_time) {
                dt = final_time - time_new;
            }

            //change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração
//#pragma unroll
            for (int i = 0; i < NEQ; i++) {
                _k_aux__[i] = _k2__[i];
                _k2__[i] = _k1__[i];
                _k1__[i] = _k_aux__[i];
            }

            //it steps the method ahead, with euler solution
//#pragma unroll
            for (int i = 0; i < NEQ; i++) {
                sv_local[i] = edos_new_euler_[i];
            }

            //verifica se o incremento para a próxima iteração ultrapassa o tempo de salvar, q neste caso é o tempo final
            if (time_new + previous_dt >= final_time) {
                //se são iguais, ja foi calculada a iteração no ultimo passo de tempo e deve-se para o laço
                //nao usar igualdade - usar esta conta, pode-se mudar a tolerância
                //printf("//d: //lf\n", threadID, fabs(final_time - time_new));
                if ((fabs(final_time - time_new) < 1.0e-5)) {
                    break;
                }
                else if (time_new < final_time) {
                    dt = previous_dt = final_time - time_new;
                    time_new += previous_dt;
                    break;
                }
                else {
                    dt = previous_dt = min_step;
                    time_new += (final_time - time_new);
                    printf("Nao era pra chegar aqui: %d: %lf\n", thread_id, final_time - time_new);
                    break;
                }
            }
            else {
                time_new +=previous_dt;
            }
        }

    }

//#pragma unroll
    for(int i = 0; i < NEQ; i++) {
        *((real*)((char*)sv + pitch * i) + thread_id) = sv_local[i];
    }

    DT = dt;
    TIME_NEW = time_new;
    PREVIOUS_DT = previous_dt;
}

#define IFNUMBER_1(name)if((celltype==1.000000000000000e+00f)) { (name) = (cmdnmax_b*1.300000000000000e+00f);    }  else{ (name) = cmdnmax_b;    }
#define IFNUMBER_2(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = 0.000000000000000e+00f;    }  else{ (name) = (5.700000000000000e-02f*expf(((-(v_old_+8.000000000000000e+01f))/6.800000000000000e+00f)));    }
#define IFNUMBER_3(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = (7.700000000000000e-01f/(1.300000000000000e-01f*(1.000000000000000e+00f+expf(((-(v_old_+1.066000000000000e+01f))/1.110000000000000e+01f)))));    }  else{ (name) = ((2.700000000000000e+00f*expf((7.900000000000000e-02f*v_old_)))+(3.100000e+05*expf((3.485000000000000e-01f*v_old_))));    }
#define IFNUMBER_4(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = 0.000000000000000e+00f;    }  else{ (name) = (((((-2.542800e+04)*expf((2.444000000000000e-01f*v_old_)))-(6.948000e-06*expf(((-4.391000000000000e-02f)*v_old_))))*(v_old_+3.778000000000000e+01f))/(1.000000000000000e+00f+expf((3.110000000000000e-01f*(v_old_+7.923000000000000e+01f)))));    }
#define IFNUMBER_5(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = ((6.000000000000000e-01f*expf((5.700000000000000e-02f*v_old_)))/(1.000000000000000e+00f+expf(((-1.000000000000000e-01f)*(v_old_+3.200000000000000e+01f)))));    }  else{ (name) = ((2.424000000000000e-02f*expf(((-1.052000000000000e-02f)*v_old_)))/(1.000000000000000e+00f+expf(((-1.378000000000000e-01f)*(v_old_+4.014000000000000e+01f)))));    }
#define IFNUMBER_6(name)if((celltype==1.000000000000000e+00f)) { (name) = (GNaL_b*6.000000000000000e-01f);    }  else{ (name) = GNaL_b;    }
#define IFNUMBER_7(name)if((celltype==1.000000000000000e+00f)) { (name) = (1.000000000000000e+00f-(9.500000000000000e-01f/(1.000000000000000e+00f+expf(((v_old_+EKshift+7.000000000000000e+01f)/5.000000000000000e+00f)))));    }  else{ (name) = 1.000000000000000e+00f;    }
#define IFNUMBER_8(name)if((celltype==1.000000000000000e+00f)) { (name) = (Gto_b*2.000000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Gto_b*2.000000000000000e+00f);    } else{ (name) = Gto_b;    }
#define IFNUMBER_9(name)if((v_old_>=3.149780000000000e+01f)) { (name) = 1.000000000000000e+00f;    }  else{ (name) = (1.076300000000000e+00f*expf(((-1.007000000000000e+00f)*expf(((-8.290000000000000e-02f)*v_old_)))));    }
#define IFNUMBER_10(name)if((celltype==1.000000000000000e+00f)) { (name) = (PCa_b*1.200000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (PCa_b*2.000000000000000e+00f);    } else{ (name) = PCa_b;    }
#define IFNUMBER_11(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKr_b*1.300000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (GKr_b*8.000000000000000e-01f);    } else{ (name) = GKr_b;    }
#define IFNUMBER_12(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKs_b*1.400000000000000e+00f);    }  else{ (name) = GKs_b;    }
#define IFNUMBER_13(name)if((celltype==1.000000000000000e+00f)) { (name) = (GK1_b*1.200000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (GK1_b*1.300000000000000e+00f);    } else{ (name) = GK1_b;    }
#define IFNUMBER_14(name)if((celltype==1.000000000000000e+00f)) { (name) = (Gncx_b*1.100000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Gncx_b*1.400000000000000e+00f);    } else{ (name) = Gncx_b;    }
#define IFNUMBER_15(name)if((celltype==1.000000000000000e+00f)) { (name) = (Pnak_b*9.000000000000000e-01f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Pnak_b*7.000000000000000e-01f);    } else{ (name) = Pnak_b;    }
#define IFNUMBER_16(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKb_b*6.000000000000000e-01f);    }  else{ (name) = GKb_b;    }
#define IFNUMBER_17(name)if((celltype==2.000000000000000e+00f)) { (name) = (calc_Jrel_inf_b*1.700000000000000e+00f);    }  else{ (name) = calc_Jrel_inf_b;    }
#define IFNUMBER_18(name)if((calc_tau_rel_b<1.000000000000000e-03f)) { (name) = 1.000000000000000e-03f;    }  else{ (name) = calc_tau_rel_b;    }
#define IFNUMBER_19(name)if((celltype==2.000000000000000e+00f)) { (name) = (calc_Jrel_infp_b*1.700000000000000e+00f);    }  else{ (name) = calc_Jrel_infp_b;    }
#define IFNUMBER_20(name)if((calc_tau_relp_b<1.000000000000000e-03f)) { (name) = 1.000000000000000e-03f;    }  else{ (name) = calc_tau_relp_b;    }
#define IFNUMBER_21(name)if((celltype==1.000000000000000e+00f)) { (name) = 1.300000000000000e+00f;    }  else{ (name) = 1.000000000000000e+00f;    }


inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt) {

    //State variables
    const real v_old_ = sv[0];
    const real CaMKt_old_ = sv[1];
    const real nai_old_ = sv[2];
    const real nass_old_ = sv[3];
    const real ki_old_ = sv[4];
    const real kss_old_ = sv[5];
    const real cai_old_ = sv[6];
    const real cass_old_ = sv[7];
    const real cansr_old_ = sv[8];
    const real cajsr_old_ = sv[9];
    const real m_old_ = sv[10];
    const real h_old_ = sv[11];
    const real j_old_ = sv[12];
    const real hp_old_ = sv[13];
    const real jp_old_ = sv[14];
    const real mL_old_ = sv[15];
    const real hL_old_ = sv[16];
    const real hLp_old_ = sv[17];
    const real a_old_ = sv[18];
    const real iF_old_ = sv[19];
    const real iS_old_ = sv[20];
    const real ap_old_ = sv[21];
    const real iFp_old_ = sv[22];
    const real iSp_old_ = sv[23];
    const real d_old_ = sv[24];
    const real ff_old_ = sv[25];
    const real fs_old_ = sv[26];
    const real fcaf_old_ = sv[27];
    const real fcas_old_ = sv[28];
    const real jca_old_ = sv[29];
    const real ffp_old_ = sv[30];
    const real fcafp_old_ = sv[31];
    const real nca_ss_old_ = sv[32];
    const real nca_i_old_ = sv[33];
    const real C3_old_ = sv[34];
    const real C2_old_ = sv[35];
    const real C1_old_ = sv[36];
    const real O_old_ = sv[37];
    const real I_old_ = sv[38];
    const real xs1_old_ = sv[39];
    const real xs2_old_ = sv[40];
    const real Jrel_np_old_ = sv[41];
    const real Jrel_p_old_ = sv[42];

//    //State variables
//    const real v_old_ =  *((real*)((char*)sv_ + pitch * 0) + threadID_);
//    const real CaMKt_old_ =  *((real*)((char*)sv_ + pitch * 1) + threadID_);
//    const real nai_old_ =  *((real*)((char*)sv_ + pitch * 2) + threadID_);
//    const real nass_old_ =  *((real*)((char*)sv_ + pitch * 3) + threadID_);
//    const real ki_old_ =  *((real*)((char*)sv_ + pitch * 4) + threadID_);
//    const real kss_old_ =  *((real*)((char*)sv_ + pitch * 5) + threadID_);
//    const real cai_old_ =  *((real*)((char*)sv_ + pitch * 6) + threadID_);
//    const real cass_old_ =  *((real*)((char*)sv_ + pitch * 7) + threadID_);
//    const real cansr_old_ =  *((real*)((char*)sv_ + pitch * 8) + threadID_);
//    const real cajsr_old_ =  *((real*)((char*)sv_ + pitch * 9) + threadID_);
//    const real m_old_ =  *((real*)((char*)sv_ + pitch * 10) + threadID_);
//    const real h_old_ =  *((real*)((char*)sv_ + pitch * 11) + threadID_);
//    const real j_old_ =  *((real*)((char*)sv_ + pitch * 12) + threadID_);
//    const real hp_old_ =  *((real*)((char*)sv_ + pitch * 13) + threadID_);
//    const real jp_old_ =  *((real*)((char*)sv_ + pitch * 14) + threadID_);
//    const real mL_old_ =  *((real*)((char*)sv_ + pitch * 15) + threadID_);
//    const real hL_old_ =  *((real*)((char*)sv_ + pitch * 16) + threadID_);
//    const real hLp_old_ =  *((real*)((char*)sv_ + pitch * 17) + threadID_);
//    const real a_old_ =  *((real*)((char*)sv_ + pitch * 18) + threadID_);
//    const real iF_old_ =  *((real*)((char*)sv_ + pitch * 19) + threadID_);
//    const real iS_old_ =  *((real*)((char*)sv_ + pitch * 20) + threadID_);
//    const real ap_old_ =  *((real*)((char*)sv_ + pitch * 21) + threadID_);
//    const real iFp_old_ =  *((real*)((char*)sv_ + pitch * 22) + threadID_);
//    const real iSp_old_ =  *((real*)((char*)sv_ + pitch * 23) + threadID_);
//    const real d_old_ =  *((real*)((char*)sv_ + pitch * 24) + threadID_);
//    const real ff_old_ =  *((real*)((char*)sv_ + pitch * 25) + threadID_);
//    const real fs_old_ =  *((real*)((char*)sv_ + pitch * 26) + threadID_);
//    const real fcaf_old_ =  *((real*)((char*)sv_ + pitch * 27) + threadID_);
//    const real fcas_old_ =  *((real*)((char*)sv_ + pitch * 28) + threadID_);
//    const real jca_old_ =  *((real*)((char*)sv_ + pitch * 29) + threadID_);
//    const real ffp_old_ =  *((real*)((char*)sv_ + pitch * 30) + threadID_);
//    const real fcafp_old_ =  *((real*)((char*)sv_ + pitch * 31) + threadID_);
//    const real nca_ss_old_ =  *((real*)((char*)sv_ + pitch * 32) + threadID_);
//    const real nca_i_old_ =  *((real*)((char*)sv_ + pitch * 33) + threadID_);
//    const real C3_old_ =  *((real*)((char*)sv_ + pitch * 34) + threadID_);
//    const real C2_old_ =  *((real*)((char*)sv_ + pitch * 35) + threadID_);
//    const real C1_old_ =  *((real*)((char*)sv_ + pitch * 36) + threadID_);
//    const real O_old_ =  *((real*)((char*)sv_ + pitch * 37) + threadID_);
//    const real I_old_ =  *((real*)((char*)sv_ + pitch * 38) + threadID_);
//    const real xs1_old_ =  *((real*)((char*)sv_ + pitch * 39) + threadID_);
//    const real xs2_old_ =  *((real*)((char*)sv_ + pitch * 40) + threadID_);
//    const real Jrel_np_old_ =  *((real*)((char*)sv_ + pitch * 41) + threadID_);
//    const real Jrel_p_old_ =  *((real*)((char*)sv_ + pitch * 42) + threadID_);

    //Parameters
    const real rad = 1.100000000000000e-03f;
    const real L = 1.000000000000000e-02f;
    const real F = 9.648500000000000e+04f;
    const real R = 8.314000000000000e+03f;
    const real T = 3.100000000000000e+02f;
    const real CaMKo = 5.000000000000000e-02f;
    const real KmCaM = 1.500000000000000e-03f;
    const real aCaMK = 5.000000000000000e-02f;
    const real bCaMK = 6.800000000000000e-04f;
    const real cmdnmax_b = 5.000000000000000e-02f;
    const real celltype = 0.000000000000000e+00f;
    const real kmcmdn = 2.380000000000000e-03f;
    const real trpnmax = 7.000000000000001e-02f;
    const real kmtrpn = 5.000000000000000e-04f;
    const real BSRmax = 4.700000000000000e-02f;
    const real KmBSR = 8.700000000000000e-04f;
    const real BSLmax = 1.124000000000000e+00f;
    const real KmBSL = 8.699999999999999e-03f;
    const real csqnmax = 1.000000000000000e+01f;
    const real kmcsqn = 8.000000000000000e-01f;
    const real zna = 1.000000000000000e+00f;
    const real nao = 1.400000000000000e+02f;
    const real zk = 1.000000000000000e+00f;
    const real ko = 5.000000000000000e+00f;
    const real PKNa = 1.833000000000000e-02f;
    const real zcl = -1.000000000000000e+00f;
    const real clo = 1.500000000000000e+02f;
    const real cli = 2.400000000000000e+01f;
    const real K_o_n = 5.000000000000000e+00f;
    const real A_atp = 2.000000000000000e+00f;
    const real K_atp = 2.500000000000000e-01f;
    const real fkatp = 0.000000000000000e+00f;
    const real gkatp = 4.319500000000000e+00f;
    const real KmCaMK = 1.500000000000000e-01f;
    const real GNa = 1.178020000000000e+01f;
    const real thL = 2.000000000000000e+02f;
    const real GNaL_b = 2.790000000000000e-02f;
    const real EKshift = 0.000000000000000e+00f;
    const real Gto_b = 1.600000000000000e-01f;
    const real offset = 0.000000000000000e+00f;
    const real vShift = 0.000000000000000e+00f;
    const real Aff = 6.000000000000000e-01f;
    const real tjca = 7.500000000000000e+01f;
    const real k2n = 5.000000000000000e+02f;
    const real Kmn = 2.000000000000000e-03f;
    const real cao = 1.800000000000000e+00f;
    const real dielConstant = 7.400000000000000e+01f;
    const real PCa_b = 8.375700000000000e-05f;
    const real ICaL_fractionSS = 8.000000000000000e-01f;
    const real beta_1 = 1.911000000000000e-01f;
    const real alpha_1 = 1.543750000000000e-01f;
    const real GKr_b = 3.210000000000000e-02f;
    const real GKs_b = 1.100000000000000e-03f;
    const real GK1_b = 6.992000000000000e-01f;
    const real qca = 1.670000000000000e-01f;
    const real qna = 5.224000000000000e-01f;
    const real kna3 = 8.812000000000000e+01f;
    const real kna1 = 1.500000000000000e+01f;
    const real kna2 = 5.000000000000000e+00f;
    const real kasymm = 1.250000000000000e+01f;
    const real kcaon = 1.500000000000000e+06f;
    const real kcaoff = 5.000000000000000e+03f;
    const real wca = 6.000000000000000e+04f;
    const real wnaca = 5.000000000000000e+03f;
    const real wna = 6.000000000000000e+04f;
    const real KmCaAct = 1.500000000000000e-04f;
    const real Gncx_b = 3.400000000000000e-03f;
    const real INaCa_fractionSS = 3.500000000000000e-01f;
    const real zca = 2.000000000000000e+00f;
    const real Knai0 = 9.073000000000000e+00f;
    const real delta = -1.550000000000000e-01f;
    const real Knao0 = 2.778000000000000e+01f;
    const real eP = 4.200000000000000e+00f;
    const real H = 1.000000000000000e-07f;
    const real Khp = 1.698000000000000e-07f;
    const real Knap = 2.240000000000000e+02f;
    const real Kxkur = 2.920000000000000e+02f;
    const real k1p = 9.495000000000000e+02f;
    const real Kki = 5.000000000000000e-01f;
    const real k1m = 1.824000000000000e+02f;
    const real MgADP = 5.000000000000000e-02f;
    const real k2p = 6.872000000000000e+02f;
    const real k2m = 3.940000000000000e+01f;
    const real Kko = 3.582000000000000e-01f;
    const real k3p = 1.899000000000000e+03f;
    const real k3m = 7.930000000000000e+04f;
    const real MgATP = 9.800000000000001e+00f;
    const real Kmgatp = 1.698000000000000e-07f;
    const real k4p = 6.390000000000000e+02f;
    const real k4m = 4.000000000000000e+01f;
    const real Pnak_b = 1.545090000000000e+01f;
    const real GKb_b = 1.890000000000000e-02f;
    const real PNab = 1.923900000000000e-09f;
    const real PCab = 5.919400000000000e-08f;
    const real GpCa = 5.000000000000000e-04f;
    const real KmCap = 5.000000000000000e-04f;
    const real Fjunc = 1.000000000000000e+00f;
    const real GClCa = 2.843000000000000e-01f;
    const real KdClCa = 1.000000000000000e-01f;
    const real GClb = 1.980000000000000e-03f;
    const real tauNa = 2.000000000000000e+00f;
    const real tauK = 2.000000000000000e+00f;
    const real tauCa = 2.000000000000000e-01f;
    const real bt = 4.750000000000000e+00f;
    const real cajsr_half = 1.700000000000000e+00f;
    const real Jrel_b = 1.537800000000000e+00f;
    const real Jup_b = 1.000000000000000e+00f;

    real calc_vcell = (1.000000000000000e+03f*3.140000000000000e+00f*rad*rad*L);	//0
    real calc_Ageo = ((2.000000000000000e+00f*3.140000000000000e+00f*rad*rad)+(2.000000000000000e+00f*3.140000000000000e+00f*rad*L));	//1

    real calc_vffrt = ((v_old_*F*F)/(R*T));	//8
    real calc_vfrt = ((v_old_*F)/(R*T));	//9
    real calc_CaMKb = ((CaMKo*(1.000000000000000e+00f-CaMKt_old_))/(1.000000000000000e+00f+(KmCaM/cass_old_)));	//11
    real calc_cmdnmax = 0.0f;
    IFNUMBER_1(calc_cmdnmax);	//14
    real calc_Bcass = (1.000000000000000e+00f/(1.000000000000000e+00f+((BSRmax*KmBSR)/powf((KmBSR+cass_old_),2.000000000000000e+00f))+((BSLmax*KmBSL)/powf((KmBSL+cass_old_),2.000000000000000e+00f))));	//21
    real calc_Bcajsr = (1.000000000000000e+00f/(1.000000000000000e+00f+((csqnmax*kmcsqn)/powf((kmcsqn+cajsr_old_),2.000000000000000e+00f))));	//24
    real calc_ENa = (((R*T)/(zna*F))*logf((nao/nai_old_)));	//26
    real calc_EK = (((R*T)/(zk*F))*logf((ko/ki_old_)));	//27
    real calc_EKs = (((R*T)/(zk*F))*logf(((ko+(PKNa*nao))/(ki_old_+(PKNa*nai_old_)))));	//28
    real calc_ECl = (((R*T)/(zcl*F))*logf((clo/cli)));	//29
    real calc_akik = powf((ko/K_o_n),2.400000000000000e-01f);	//30
    real calc_bkik = (1.000000000000000e+00f/(1.000000000000000e+00f+powf((A_atp/K_atp),2.000000000000000e+00f)));	//31
    real calc_mss = (1.000000000000000e+00f/powf((1.000000000000000e+00f+expf(((-(v_old_+5.686000000000000e+01f))/9.029999999999999e+00f))),2.000000000000000e+00f));	//33
    real calc_tm = ((1.292000000000000e-01f*expf((-powf(((v_old_+4.579000000000000e+01f)/1.554000000000000e+01f),2.000000000000000e+00f))))+(6.487000000000000e-02f*expf((-powf(((v_old_-4.823000000000000e+00f)/5.112000000000000e+01f),2.000000000000000e+00f)))));	//34
    real calc_hss = (1.000000000000000e+00f/powf((1.000000000000000e+00f+expf(((v_old_+7.155000000000000e+01f)/7.430000000000000e+00f))),2.000000000000000e+00f));	//36
    real calc_ah = 0.0f;
    IFNUMBER_2(calc_ah);	//37
    real calc_bh = 0.0f;
    IFNUMBER_3(calc_bh);	//38
    real calc_aj = 0.0f;
    IFNUMBER_4(calc_aj);	//41
    real calc_bj = 0.0f;
    IFNUMBER_5(calc_bj);	//42
    real calc_hssp = (1.000000000000000e+00f/powf((1.000000000000000e+00f+expf(((v_old_+7.755000000000000e+01f)/7.430000000000000e+00f))),2.000000000000000e+00f));	//46
    real calc_mLss = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((-(v_old_+4.285000000000000e+01f))/5.264000000000000e+00f))));	//52
    real calc_tmL = ((1.292000000000000e-01f*expf((-powf(((v_old_+4.579000000000000e+01f)/1.554000000000000e+01f),2.000000000000000e+00f))))+(6.487000000000000e-02f*expf((-powf(((v_old_-4.823000000000000e+00f)/5.112000000000000e+01f),2.000000000000000e+00f)))));	//53
    real calc_hLss = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((v_old_+8.761000000000000e+01f)/7.488000000000000e+00f))));	//55
    real calc_hLssp = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((v_old_+9.381000000000000e+01f)/7.488000000000000e+00f))));	//57
    real calc_thLp = (3.000000000000000e+00f*thL);	//58
    real calc_GNaL = 0.0f;
    IFNUMBER_6(calc_GNaL);	//60
    real calc_ass = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((-((v_old_+EKshift)-1.434000000000000e+01f))/1.482000000000000e+01f))));	//63
    real calc_ta = (1.051500000000000e+00f/((1.000000000000000e+00f/(1.208900000000000e+00f*(1.000000000000000e+00f+expf(((-((v_old_+EKshift)-1.840990000000000e+01f))/2.938140000000000e+01f)))))+(3.500000000000000e+00f/(1.000000000000000e+00f+expf(((v_old_+EKshift+1.000000000000000e+02f)/2.938140000000000e+01f))))));	//64
    real calc_iss = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((v_old_+EKshift+4.394000000000000e+01f)/5.711000000000000e+00f))));	//66
    real calc_delta_epi = 0.0f;
    IFNUMBER_7(calc_delta_epi);	//67
    real calc_tiF_b = (4.562000000000000e+00f+(1.000000000000000e+00f/((3.933000000000000e-01f*expf(((-(v_old_+EKshift+1.000000000000000e+02f))/1.000000000000000e+02f)))+(8.004000000000000e-02f*expf(((v_old_+EKshift+5.000000000000000e+01f)/1.659000000000000e+01f))))));	//68
    real calc_tiS_b = (2.362000000000000e+01f+(1.000000000000000e+00f/((1.416000000000000e-03f*expf(((-(v_old_+EKshift+9.652000000000000e+01f))/5.905000000000000e+01f)))+(1.780000e-08*expf(((v_old_+EKshift+1.141000000000000e+02f)/8.079000000000001e+00f))))));	//69
    real calc_AiF = (1.000000000000000e+00f/(1.000000000000000e+00f+expf((((v_old_+EKshift)-2.136000000000000e+02f)/1.512000000000000e+02f))));	//72
    real calc_assp = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((-((v_old_+EKshift)-2.434000000000000e+01f))/1.482000000000000e+01f))));	//77
    real calc_dti_develop = (1.354000000000000e+00f+(1.000000e-04/(expf((((v_old_+EKshift)-1.674000000000000e+02f)/1.589000000000000e+01f))+expf(((-((v_old_+EKshift)-1.223000000000000e+01f))/2.154000000000000e-01f)))));	//79
    real calc_dti_recover = (1.000000000000000e+00f-(5.000000000000000e-01f/(1.000000000000000e+00f+expf(((v_old_+EKshift+7.000000000000000e+01f)/2.000000000000000e+01f)))));	//80
    real calc_Gto = 0.0f;
    IFNUMBER_8(calc_Gto);	//86
    real calc_dss = 0.0f;
    IFNUMBER_9(calc_dss);	//89
    real calc_td = (offset+6.000000000000000e-01f+(1.000000000000000e+00f/(expf(((-5.000000000000000e-02f)*(v_old_+vShift+6.000000000000000e+00f)))+expf((9.000000000000000e-02f*(v_old_+vShift+1.400000000000000e+01f))))));	//90
    real calc_fss = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((v_old_+1.958000000000000e+01f)/3.696000000000000e+00f))));	//92
    real calc_tff = (7.000000000000000e+00f+(1.000000000000000e+00f/((4.500000000000000e-03f*expf(((-(v_old_+2.000000000000000e+01f))/1.000000000000000e+01f)))+(4.500000000000000e-03f*expf(((v_old_+2.000000000000000e+01f)/1.000000000000000e+01f))))));	//93
    real calc_tfs = (1.000000000000000e+03f+(1.000000000000000e+00f/((3.500000000000000e-05f*expf(((-(v_old_+5.000000000000000e+00f))/4.000000000000000e+00f)))+(3.500000000000000e-05f*expf(((v_old_+5.000000000000000e+00f)/6.000000000000000e+00f))))));	//94
    real calc_Afs = (1.000000000000000e+00f-Aff);	//95
    real calc_tfcaf = (7.000000000000000e+00f+(1.000000000000000e+00f/((4.000000000000000e-02f*expf(((-(v_old_-4.000000000000000e+00f))/7.000000000000000e+00f)))+(4.000000000000000e-02f*expf(((v_old_-4.000000000000000e+00f)/7.000000000000000e+00f))))));	//100
    real calc_tfcas = (1.000000000000000e+02f+(1.000000000000000e+00f/((1.200000000000000e-04f*expf(((-v_old_)/3.000000000000000e+00f)))+(1.200000000000000e-04f*expf((v_old_/7.000000000000000e+00f))))));	//101
    real calc_Afcaf = (3.000000000000000e-01f+(6.000000000000000e-01f/(1.000000000000000e+00f+expf(((v_old_-1.000000000000000e+01f)/1.000000000000000e+01f)))));	//102
    real calc_jcass = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((v_old_+1.808000000000000e+01f)/2.791600000000000e+00f))));	//107
    real calc_km2n = (jca_old_*1.000000000000000e+00f);	//115
    real calc_Io = ((5.000000000000000e-01f*(nao+ko+clo+(4.000000000000000e+00f*cao)))/1.000000000000000e+03f);	//118
    real calc_Iss = ((5.000000000000000e-01f*(nass_old_+kss_old_+cli+(4.000000000000000e+00f*cass_old_)))/1.000000000000000e+03f);	//119
    real calc_constA = (1.820000e+06*powf((dielConstant*T),(-1.500000000000000e+00f)));	//120
    real calc_PCa = 0.0f;
    IFNUMBER_10(calc_PCa);	//130
    real calc_Ii = ((5.000000000000000e-01f*(nai_old_+ki_old_+cli+(4.000000000000000e+00f*cai_old_)))/1.000000000000000e+03f);	//142
    real calc_GKr = 0.0f;
    IFNUMBER_11(calc_GKr);	//168
    real calc_xs1ss = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((-(v_old_+1.160000000000000e+01f))/8.932000000000000e+00f))));	//170
    real calc_txs1 = (8.173000000000000e+02f+(1.000000000000000e+00f/((2.326000e-04*expf(((v_old_+4.828000000000000e+01f)/1.780000000000000e+01f)))+(1.292000000000000e-03f*expf(((-(v_old_+2.100000000000000e+02f))/2.300000000000000e+02f))))));	//171
    real calc_txs2 = (1.000000000000000e+00f/((1.000000000000000e-02f*expf(((v_old_-5.000000000000000e+01f)/2.000000000000000e+01f)))+(1.930000000000000e-02f*expf(((-(v_old_+6.654000000000001e+01f))/3.100000000000000e+01f)))));	//174
    real calc_KsCa = (1.000000000000000e+00f+(6.000000000000000e-01f/(1.000000000000000e+00f+powf((3.800000e-05/cai_old_),1.400000000000000e+00f))));	//176
    real calc_GKs = 0.0f;
    IFNUMBER_12(calc_GKs);	//177
    real calc_GK1 = 0.0f;
    IFNUMBER_13(calc_GK1);	//182
    real calc_h4_i = (1.000000000000000e+00f+((nai_old_/kna1)*(1.000000000000000e+00f+(nai_old_/kna2))));	//189
    real calc_h10_i = (kasymm+1.000000000000000e+00f+((nao/kna1)*(1.000000000000000e+00f+(nao/kna2))));	//195
    real calc_k2_i = kcaoff;	//199
    real calc_k5_i = kcaoff;	//206
    real calc_allo_i = (1.000000000000000e+00f/(1.000000000000000e+00f+powf((KmCaAct/cai_old_),2.000000000000000e+00f)));	//218
    real calc_Gncx = 0.0f;
    IFNUMBER_14(calc_Gncx);	//221
    real calc_h4_ss = (1.000000000000000e+00f+((nass_old_/kna1)*(1.000000000000000e+00f+(nass_old_/kna2))));	//226
    real calc_h10_ss = (kasymm+1.000000000000000e+00f+((nao/kna1)*(1.000000000000000e+00f+(nao/kna2))));	//232
    real calc_k2_ss = kcaoff;	//236
    real calc_k5_ss = kcaoff;	//243
    real calc_allo_ss = (1.000000000000000e+00f/(1.000000000000000e+00f+powf((KmCaAct/cass_old_),2.000000000000000e+00f)));	//255
    real calc_P = (eP/(1.000000000000000e+00f+(H/Khp)+(nai_old_/Knap)+(ki_old_/Kxkur)));	//261
    real calc_b1 = (k1m*MgADP);	//263
    real calc_a2 = k2p;	//264
    real calc_a4 = (((k4p*MgATP)/Kmgatp)/(1.000000000000000e+00f+(MgATP/Kmgatp)));	//268
    real calc_Pnak = 0.0f;
    IFNUMBER_15(calc_Pnak);	//280
    real calc_xkb = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((-(v_old_-1.089680000000000e+01f))/2.398710000000000e+01f))));	//282
    real calc_GKb = 0.0f;
    IFNUMBER_16(calc_GKb);	//283
    real calc_IpCa = ((GpCa*cai_old_)/(KmCap+cai_old_));	//287
    real calc_JdiffNa = ((nass_old_-nai_old_)/tauNa);	//292
    real calc_JdiffK = ((kss_old_-ki_old_)/tauK);	//293
    real calc_Jdiff = ((cass_old_-cai_old_)/tauCa);	//294
    real calc_a_rel = ((5.000000000000000e-01f*bt)/1.000000000000000e+00f);	//295
    real calc_tau_rel_b = (bt/(1.000000000000000e+00f+(1.230000000000000e-02f/cajsr_old_)));	//298
    real calc_btp = (1.250000000000000e+00f*bt);	//301
    real calc_upScale = 0.0f;
    IFNUMBER_21(calc_upScale);	//310
    real calc_Jleak = ((4.882500000000000e-03f*cansr_old_)/1.500000000000000e+01f);	//314
    real calc_Jtr = ((cansr_old_-cajsr_old_)/6.000000000000000e+01f);	//316
    real calc_Acap = (2.000000000000000e+00f*calc_Ageo);	//2
    real calc_vmyo = (6.800000000000000e-01f*calc_vcell);	//3
    real calc_vnsr = (5.520000000000000e-02f*calc_vcell);	//4
    real calc_vjsr = (4.800000000000000e-03f*calc_vcell);	//5
    real calc_vss = (2.000000000000000e-02f*calc_vcell);	//6
    real calc_CaMKa = (calc_CaMKb+CaMKt_old_);	//12
    real calc_Bcai = (1.000000000000000e+00f/(1.000000000000000e+00f+((calc_cmdnmax*kmcmdn)/powf((kmcmdn+cai_old_),2.000000000000000e+00f))+((trpnmax*kmtrpn)/powf((kmtrpn+cai_old_),2.000000000000000e+00f))));	//19
    real calc_I_katp = (fkatp*gkatp*calc_akik*calc_bkik*(v_old_-calc_EK));	//32
    real calc_th = (1.000000000000000e+00f/(calc_ah+calc_bh));	//39
    real calc_jss = calc_hss;	//43
    real calc_tj = (1.000000000000000e+00f/(calc_aj+calc_bj));	//44
    real calc_tiF = (calc_tiF_b*calc_delta_epi);	//70
    real calc_tiS = (calc_tiS_b*calc_delta_epi);	//71
    real calc_AiS = (1.000000000000000e+00f-calc_AiF);	//73
    real calc_f = ((Aff*ff_old_)+(calc_Afs*fs_old_));	//98
    real calc_fcass = calc_fss;	//99
    real calc_Afcas = (1.000000000000000e+00f-calc_Afcaf);	//103
    real calc_tffp = (2.500000000000000e+00f*calc_tff);	//109
    real calc_fp = ((Aff*ffp_old_)+(calc_Afs*fs_old_));	//111
    real calc_tfcafp = (2.500000000000000e+00f*calc_tfcaf);	//112
    real calc_anca_ss = (1.000000000000000e+00f/((k2n/calc_km2n)+powf((1.000000000000000e+00f+(Kmn/cass_old_)),4.000000000000000e+00f)));	//116
    real calc_gamma_cass = expf(((-calc_constA)*4.000000000000000e+00f*((powf(calc_Iss,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Iss,1.0/2.0)))-(3.000000000000000e-01f*calc_Iss))));	//121
    real calc_gamma_cao = expf(((-calc_constA)*4.000000000000000e+00f*((powf(calc_Io,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Io,1.0/2.0)))-(3.000000000000000e-01f*calc_Io))));	//122
    real calc_gamma_nass = expf(((-calc_constA)*1.000000000000000e+00f*((powf(calc_Iss,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Iss,1.0/2.0)))-(3.000000000000000e-01f*calc_Iss))));	//123
    real calc_gamma_nao = expf(((-calc_constA)*1.000000000000000e+00f*((powf(calc_Io,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Io,1.0/2.0)))-(3.000000000000000e-01f*calc_Io))));	//124
    real calc_gamma_kss = expf(((-calc_constA)*1.000000000000000e+00f*((powf(calc_Iss,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Iss,1.0/2.0)))-(3.000000000000000e-01f*calc_Iss))));	//125
    real calc_gamma_ko = expf(((-calc_constA)*1.000000000000000e+00f*((powf(calc_Io,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Io,1.0/2.0)))-(3.000000000000000e-01f*calc_Io))));	//126
    real calc_PCap = (1.100000000000000e+00f*calc_PCa);	//131
    real calc_PCaNa = (1.250000000000000e-03f*calc_PCa);	//132
    real calc_PCaK = (3.574000e-04*calc_PCa);	//133
    real calc_anca_i = (1.000000000000000e+00f/((k2n/calc_km2n)+powf((1.000000000000000e+00f+(Kmn/cai_old_)),4.000000000000000e+00f)));	//140
    real calc_gamma_cai = expf(((-calc_constA)*4.000000000000000e+00f*((powf(calc_Ii,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Ii,1.0/2.0)))-(3.000000000000000e-01f*calc_Ii))));	//143
    real calc_gamma_nai = expf(((-calc_constA)*1.000000000000000e+00f*((powf(calc_Ii,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Ii,1.0/2.0)))-(3.000000000000000e-01f*calc_Ii))));	//144
    real calc_gamma_ki = expf(((-calc_constA)*1.000000000000000e+00f*((powf(calc_Ii,1.0/2.0)/(1.000000000000000e+00f+powf(calc_Ii,1.0/2.0)))-(3.000000000000000e-01f*calc_Ii))));	//145
    real calc_alpha = (1.161000000000000e-01f*expf((2.990000000000000e-01f*calc_vfrt)));	//155
    real calc_beta = (2.442000000000000e-01f*expf(((-1.604000000000000e+00f)*calc_vfrt)));	//156
    real calc_alpha_2 = (5.780000000000000e-02f*expf((9.710000000000000e-01f*calc_vfrt)));	//157
    real calc_beta_2 = (0.349000e-03*expf(((-1.062000000000000e+00f)*calc_vfrt)));	//158
    real calc_alpha_i = (2.533000000000000e-01f*expf((5.953000000000001e-01f*calc_vfrt)));	//159
    real calc_beta_i = (6.525000000000000e-02f*expf(((-8.209000000000000e-01f)*calc_vfrt)));	//160
    real calc_alpha_C2ToI = (0.520000e-04*expf((1.525000000000000e+00f*calc_vfrt)));	//161
    real calc_IKr = (calc_GKr*powf((ko/5.000000000000000e+00f),1.0/2.0)*O_old_*(v_old_-calc_EK));	//169
    real calc_xs2ss = calc_xs1ss;	//173
    real calc_IKs = (calc_GKs*calc_KsCa*xs1_old_*xs2_old_*(v_old_-calc_EKs));	//178
    real calc_aK1 = (4.094000000000000e+00f/(1.000000000000000e+00f+expf((1.217000000000000e-01f*((v_old_-calc_EK)-4.993400000000000e+01f)))));	//179
    real calc_bK1 = (((1.572000000000000e+01f*expf((6.740000000000000e-02f*((v_old_-calc_EK)-3.257000000000000e+00f))))+expf((6.180000000000000e-02f*((v_old_-calc_EK)-5.943099999999999e+02f))))/(1.000000000000000e+00f+expf(((-1.629000000000000e-01f)*((v_old_-calc_EK)+1.420700000000000e+01f)))));	//180
    real calc_hca = expf((qca*calc_vfrt));	//184
    real calc_hna = expf((qna*calc_vfrt));	//185
    real calc_h5_i = ((nai_old_*nai_old_)/(calc_h4_i*kna1*kna2));	//190
    real calc_h6_i = (1.000000000000000e+00f/calc_h4_i);	//191
    real calc_h11_i = ((nao*nao)/(calc_h10_i*kna1*kna2));	//196
    real calc_h12_i = (1.000000000000000e+00f/calc_h10_i);	//197
    real calc_h5_ss = ((nass_old_*nass_old_)/(calc_h4_ss*kna1*kna2));	//227
    real calc_h6_ss = (1.000000000000000e+00f/calc_h4_ss);	//228
    real calc_h11_ss = ((nao*nao)/(calc_h10_ss*kna1*kna2));	//233
    real calc_h12_ss = (1.000000000000000e+00f/calc_h10_ss);	//234
    real calc_Knai = (Knai0*expf(((delta*calc_vfrt)/3.000000000000000e+00f)));	//259
    real calc_Knao = (Knao0*expf((((1.000000000000000e+00f-delta)*calc_vfrt)/3.000000000000000e+00f)));	//260
    real calc_b3 = ((k3m*calc_P*H)/(1.000000000000000e+00f+(MgATP/Kmgatp)));	//267
    real calc_IKb = (calc_GKb*calc_xkb*(v_old_-calc_EK));	//284
    real calc_INab = ((PNab*calc_vffrt*((nai_old_*expf(calc_vfrt))-nao))/(expf(calc_vfrt)-1.000000000000000e+00f));	//285
    real calc_IClCa_junc = (((Fjunc*GClCa)/(1.000000000000000e+00f+(KdClCa/cass_old_)))*(v_old_-calc_ECl));	//288
    real calc_IClCa_sl = ((((1.000000000000000e+00f-Fjunc)*GClCa)/(1.000000000000000e+00f+(KdClCa/cai_old_)))*(v_old_-calc_ECl));	//289
    real calc_IClb = (GClb*(v_old_-calc_ECl));	//291
    real calc_tau_rel = 0.0f;
    IFNUMBER_18(calc_tau_rel);	//299
    real calc_a_relp = ((5.000000000000000e-01f*calc_btp)/1.000000000000000e+00f);	//302
    real calc_tau_relp_b = (calc_btp/(1.000000000000000e+00f+(1.230000000000000e-02f/cajsr_old_)));	//305
    real calc_Jupnp = ((calc_upScale*5.425000000000000e-03f*cai_old_)/(cai_old_+9.200000000000000e-04f));	//311
    real calc_Jupp = ((calc_upScale*2.750000000000000e+00f*5.425000000000000e-03f*cai_old_)/((cai_old_+9.200000000000000e-04f)-1.700000000000000e-04f));	//312
    real calc_tjp = (1.460000000000000e+00f*calc_tj);	//48
    real calc_fINap = (1.000000000000000e+00f/(1.000000000000000e+00f+(KmCaMK/calc_CaMKa)));	//50
    real calc_fINaLp = (1.000000000000000e+00f/(1.000000000000000e+00f+(KmCaMK/calc_CaMKa)));	//61
    real calc_i = ((calc_AiF*iF_old_)+(calc_AiS*iS_old_));	//76
    real calc_tiFp = (calc_dti_develop*calc_dti_recover*calc_tiF);	//81
    real calc_tiSp = (calc_dti_develop*calc_dti_recover*calc_tiS);	//82
    real calc_ip = ((calc_AiF*iFp_old_)+(calc_AiS*iSp_old_));	//85
    real calc_fItop = (1.000000000000000e+00f/(1.000000000000000e+00f+(KmCaMK/calc_CaMKa)));	//87
    real calc_fca = ((calc_Afcaf*fcaf_old_)+(calc_Afcas*fcas_old_));	//106
    real calc_fcap = ((calc_Afcaf*fcafp_old_)+(calc_Afcas*fcas_old_));	//114
    real calc_PhiCaL_ss = ((4.000000000000000e+00f*calc_vffrt*((calc_gamma_cass*cass_old_*expf((2.000000000000000e+00f*calc_vfrt)))-(calc_gamma_cao*cao)))/(expf((2.000000000000000e+00f*calc_vfrt))-1.000000000000000e+00f));	//127
    real calc_PhiCaNa_ss = ((1.000000000000000e+00f*calc_vffrt*((calc_gamma_nass*nass_old_*expf((1.000000000000000e+00f*calc_vfrt)))-(calc_gamma_nao*nao)))/(expf((1.000000000000000e+00f*calc_vfrt))-1.000000000000000e+00f));	//128
    real calc_PhiCaK_ss = ((1.000000000000000e+00f*calc_vffrt*((calc_gamma_kss*kss_old_*expf((1.000000000000000e+00f*calc_vfrt)))-(calc_gamma_ko*ko)))/(expf((1.000000000000000e+00f*calc_vfrt))-1.000000000000000e+00f));	//129
    real calc_PCaNap = (1.250000000000000e-03f*calc_PCap);	//134
    real calc_PCaKp = (3.574000e-04*calc_PCap);	//135
    real calc_fICaLp = (1.000000000000000e+00f/(1.000000000000000e+00f+(KmCaMK/calc_CaMKa)));	//136
    real calc_PhiCaL_i = ((4.000000000000000e+00f*calc_vffrt*((calc_gamma_cai*cai_old_*expf((2.000000000000000e+00f*calc_vfrt)))-(calc_gamma_cao*cao)))/(expf((2.000000000000000e+00f*calc_vfrt))-1.000000000000000e+00f));	//146
    real calc_PhiCaNa_i = ((1.000000000000000e+00f*calc_vffrt*((calc_gamma_nai*nai_old_*expf((1.000000000000000e+00f*calc_vfrt)))-(calc_gamma_nao*nao)))/(expf((1.000000000000000e+00f*calc_vfrt))-1.000000000000000e+00f));	//147
    real calc_PhiCaK_i = ((1.000000000000000e+00f*calc_vffrt*((calc_gamma_ki*ki_old_*expf((1.000000000000000e+00f*calc_vfrt)))-(calc_gamma_ko*ko)))/(expf((1.000000000000000e+00f*calc_vfrt))-1.000000000000000e+00f));	//148
    real calc_beta_ItoC2 = ((calc_beta_2*calc_beta_i*calc_alpha_C2ToI)/(calc_alpha_2*calc_alpha_i));	//162
    real calc_K1ss = (calc_aK1/(calc_aK1+calc_bK1));	//181
    real calc_h1_i = (1.000000000000000e+00f+((nai_old_/kna3)*(1.000000000000000e+00f+calc_hna)));	//186
    real calc_h7_i = (1.000000000000000e+00f+((nao/kna3)*(1.000000000000000e+00f+(1.000000000000000e+00f/calc_hna))));	//192
    real calc_k1_i = (calc_h12_i*cao*kcaon);	//198
    real calc_k6_i = (calc_h6_i*cai_old_*kcaon);	//207
    real calc_h1_ss = (1.000000000000000e+00f+((nass_old_/kna3)*(1.000000000000000e+00f+calc_hna)));	//223
    real calc_h7_ss = (1.000000000000000e+00f+((nao/kna3)*(1.000000000000000e+00f+(1.000000000000000e+00f/calc_hna))));	//229
    real calc_k1_ss = (calc_h12_ss*cao*kcaon);	//235
    real calc_k6_ss = (calc_h6_ss*cass_old_*kcaon);	//244
    real calc_a1 = ((k1p*powf((nai_old_/calc_Knai),3.000000000000000e+00f))/((powf((1.000000000000000e+00f+(nai_old_/calc_Knai)),3.000000000000000e+00f)+powf((1.000000000000000e+00f+(ki_old_/Kki)),2.000000000000000e+00f))-1.000000000000000e+00f));	//262
    real calc_b2 = ((k2m*powf((nao/calc_Knao),3.000000000000000e+00f))/((powf((1.000000000000000e+00f+(nao/calc_Knao)),3.000000000000000e+00f)+powf((1.000000000000000e+00f+(ko/Kko)),2.000000000000000e+00f))-1.000000000000000e+00f));	//265
    real calc_a3 = ((k3p*powf((ko/Kko),2.000000000000000e+00f))/((powf((1.000000000000000e+00f+(nao/calc_Knao)),3.000000000000000e+00f)+powf((1.000000000000000e+00f+(ko/Kko)),2.000000000000000e+00f))-1.000000000000000e+00f));	//266
    real calc_b4 = ((k4m*powf((ki_old_/Kki),2.000000000000000e+00f))/((powf((1.000000000000000e+00f+(nai_old_/calc_Knai)),3.000000000000000e+00f)+powf((1.000000000000000e+00f+(ki_old_/Kki)),2.000000000000000e+00f))-1.000000000000000e+00f));	//269
    real calc_ICab = ((PCab*4.000000000000000e+00f*calc_vffrt*((calc_gamma_cai*cai_old_*expf((2.000000000000000e+00f*calc_vfrt)))-(calc_gamma_cao*cao)))/(expf((2.000000000000000e+00f*calc_vfrt))-1.000000000000000e+00f));	//286
    real calc_IClCa = (calc_IClCa_junc+calc_IClCa_sl);	//290
    real calc_tau_relp = 0.0f;
    IFNUMBER_20(calc_tau_relp);	//306
    real calc_fJrelp = (1.000000000000000e+00f/(1.000000000000000e+00f+(KmCaMK/calc_CaMKa)));	//308
    real calc_fJupp = (1.000000000000000e+00f/(1.000000000000000e+00f+(KmCaMK/calc_CaMKa)));	//313
    real calc_INa = (GNa*(v_old_-calc_ENa)*powf(m_old_,3.000000000000000e+00f)*(((1.000000000000000e+00f-calc_fINap)*h_old_*j_old_)+(calc_fINap*hp_old_*jp_old_)));	//51
    real calc_INaL = (calc_GNaL*(v_old_-calc_ENa)*mL_old_*(((1.000000000000000e+00f-calc_fINaLp)*hL_old_)+(calc_fINaLp*hLp_old_)));	//62
    real calc_Ito = (calc_Gto*(v_old_-calc_EK)*(((1.000000000000000e+00f-calc_fItop)*a_old_*calc_i)+(calc_fItop*ap_old_*calc_ip)));	//88
    real calc_ICaL_ss = (ICaL_fractionSS*(((1.000000000000000e+00f-calc_fICaLp)*calc_PCa*calc_PhiCaL_ss*d_old_*((calc_f*(1.000000000000000e+00f-nca_ss_old_))+(jca_old_*calc_fca*nca_ss_old_)))+(calc_fICaLp*calc_PCap*calc_PhiCaL_ss*d_old_*((calc_fp*(1.000000000000000e+00f-nca_ss_old_))+(jca_old_*calc_fcap*nca_ss_old_)))));	//137
    real calc_ICaNa_ss = (ICaL_fractionSS*(((1.000000000000000e+00f-calc_fICaLp)*calc_PCaNa*calc_PhiCaNa_ss*d_old_*((calc_f*(1.000000000000000e+00f-nca_ss_old_))+(jca_old_*calc_fca*nca_ss_old_)))+(calc_fICaLp*calc_PCaNap*calc_PhiCaNa_ss*d_old_*((calc_fp*(1.000000000000000e+00f-nca_ss_old_))+(jca_old_*calc_fcap*nca_ss_old_)))));	//138
    real calc_ICaK_ss = (ICaL_fractionSS*(((1.000000000000000e+00f-calc_fICaLp)*calc_PCaK*calc_PhiCaK_ss*d_old_*((calc_f*(1.000000000000000e+00f-nca_ss_old_))+(jca_old_*calc_fca*nca_ss_old_)))+(calc_fICaLp*calc_PCaKp*calc_PhiCaK_ss*d_old_*((calc_fp*(1.000000000000000e+00f-nca_ss_old_))+(jca_old_*calc_fcap*nca_ss_old_)))));	//139
    real calc_ICaL_i = ((1.000000000000000e+00f-ICaL_fractionSS)*(((1.000000000000000e+00f-calc_fICaLp)*calc_PCa*calc_PhiCaL_i*d_old_*((calc_f*(1.000000000000000e+00f-nca_i_old_))+(jca_old_*calc_fca*nca_i_old_)))+(calc_fICaLp*calc_PCap*calc_PhiCaL_i*d_old_*((calc_fp*(1.000000000000000e+00f-nca_i_old_))+(jca_old_*calc_fcap*nca_i_old_)))));	//149
    real calc_ICaNa_i = ((1.000000000000000e+00f-ICaL_fractionSS)*(((1.000000000000000e+00f-calc_fICaLp)*calc_PCaNa*calc_PhiCaNa_i*d_old_*((calc_f*(1.000000000000000e+00f-nca_i_old_))+(jca_old_*calc_fca*nca_i_old_)))+(calc_fICaLp*calc_PCaNap*calc_PhiCaNa_i*d_old_*((calc_fp*(1.000000000000000e+00f-nca_i_old_))+(jca_old_*calc_fcap*nca_i_old_)))));	//150
    real calc_ICaK_i = ((1.000000000000000e+00f-ICaL_fractionSS)*(((1.000000000000000e+00f-calc_fICaLp)*calc_PCaK*calc_PhiCaK_i*d_old_*((calc_f*(1.000000000000000e+00f-nca_i_old_))+(jca_old_*calc_fca*nca_i_old_)))+(calc_fICaLp*calc_PCaKp*calc_PhiCaK_i*d_old_*((calc_fp*(1.000000000000000e+00f-nca_i_old_))+(jca_old_*calc_fcap*nca_i_old_)))));	//151
    real calc_IK1 = (calc_GK1*powf((ko/5.000000000000000e+00f),1.0/2.0)*calc_K1ss*(v_old_-calc_EK));	//183
    real calc_h2_i = ((nai_old_*calc_hna)/(kna3*calc_h1_i));	//187
    real calc_h3_i = (1.000000000000000e+00f/calc_h1_i);	//188
    real calc_h8_i = (nao/(kna3*calc_hna*calc_h7_i));	//193
    real calc_h9_i = (1.000000000000000e+00f/calc_h7_i);	//194
    real calc_h2_ss = ((nass_old_*calc_hna)/(kna3*calc_h1_ss));	//224
    real calc_h3_ss = (1.000000000000000e+00f/calc_h1_ss);	//225
    real calc_h8_ss = (nao/(kna3*calc_hna*calc_h7_ss));	//230
    real calc_h9_ss = (1.000000000000000e+00f/calc_h7_ss);	//231
    real calc_x1 = ((calc_a4*calc_a1*calc_a2)+(calc_b2*calc_b4*calc_b3)+(calc_a2*calc_b4*calc_b3)+(calc_b3*calc_a1*calc_a2));	//270
    real calc_x2 = ((calc_b2*calc_b1*calc_b4)+(calc_a1*calc_a2*calc_a3)+(calc_a3*calc_b1*calc_b4)+(calc_a2*calc_a3*calc_b4));	//271
    real calc_x3 = ((calc_a2*calc_a3*calc_a4)+(calc_b3*calc_b2*calc_b1)+(calc_b2*calc_b1*calc_a4)+(calc_a3*calc_a4*calc_b1));	//272
    real calc_x4 = ((calc_b4*calc_b3*calc_b2)+(calc_a3*calc_a4*calc_a1)+(calc_b2*calc_a4*calc_a1)+(calc_b3*calc_b2*calc_a1));	//273
    real calc_Jrel = (Jrel_b*(((1.000000000000000e+00f-calc_fJrelp)*Jrel_np_old_)+(calc_fJrelp*Jrel_p_old_)));	//309
    real calc_Jup = (Jup_b*((((1.000000000000000e+00f-calc_fJupp)*calc_Jupnp)+(calc_fJupp*calc_Jupp))-calc_Jleak));	//315
    real calc_ICaL = (calc_ICaL_ss+calc_ICaL_i);	//152
    real calc_ICaNa = (calc_ICaNa_ss+calc_ICaNa_i);	//153
    real calc_ICaK = (calc_ICaK_ss+calc_ICaK_i);	//154
    real calc_k3p_i = (calc_h9_i*wca);	//200
    real calc_k3pp_i = (calc_h8_i*wnaca);	//201
    real calc_k4p_i = ((calc_h3_i*wca)/calc_hca);	//203
    real calc_k4pp_i = (calc_h2_i*wnaca);	//204
    real calc_k7_i = (calc_h5_i*calc_h2_i*wna);	//208
    real calc_k8_i = (calc_h8_i*calc_h11_i*wna);	//209
    real calc_k3p_ss = (calc_h9_ss*wca);	//237
    real calc_k3pp_ss = (calc_h8_ss*wnaca);	//238
    real calc_k4p_ss = ((calc_h3_ss*wca)/calc_hca);	//240
    real calc_k4pp_ss = (calc_h2_ss*wnaca);	//241
    real calc_k7_ss = (calc_h5_ss*calc_h2_ss*wna);	//245
    real calc_k8_ss = (calc_h8_ss*calc_h11_ss*wna);	//246
    real calc_E1 = (calc_x1/(calc_x1+calc_x2+calc_x3+calc_x4));	//274
    real calc_E2 = (calc_x2/(calc_x1+calc_x2+calc_x3+calc_x4));	//275
    real calc_E3 = (calc_x3/(calc_x1+calc_x2+calc_x3+calc_x4));	//276
    real calc_E4 = (calc_x4/(calc_x1+calc_x2+calc_x3+calc_x4));	//277
    real calc_Jrel_inf_b = ((((-calc_a_rel)*calc_ICaL_ss)/1.000000000000000e+00f)/(1.000000000000000e+00f+powf((cajsr_half/cajsr_old_),8.000000000000000e+00f)));	//296
    real calc_Jrel_infp_b = ((((-calc_a_relp)*calc_ICaL_ss)/1.000000000000000e+00f)/(1.000000000000000e+00f+powf((cajsr_half/cajsr_old_),8.000000000000000e+00f)));	//303
    real calc_k3_i = (calc_k3p_i+calc_k3pp_i);	//202
    real calc_k4_i = (calc_k4p_i+calc_k4pp_i);	//205
    real calc_k3_ss = (calc_k3p_ss+calc_k3pp_ss);	//239
    real calc_k4_ss = (calc_k4p_ss+calc_k4pp_ss);	//242
    real calc_JnakNa = (3.000000000000000e+00f*((calc_E1*calc_a3)-(calc_E2*calc_b3)));	//278
    real calc_JnakK = (2.000000000000000e+00f*((calc_E4*calc_b1)-(calc_E3*calc_a1)));	//279
    real calc_Jrel_inf = 0.0f;
    IFNUMBER_17(calc_Jrel_inf);	//297
    real calc_Jrel_infp = 0.0f;
    IFNUMBER_19(calc_Jrel_infp);	//304
    real calc_x1_i = ((calc_k2_i*calc_k4_i*(calc_k7_i+calc_k6_i))+(calc_k5_i*calc_k7_i*(calc_k2_i+calc_k3_i)));	//210
    real calc_x2_i = ((calc_k1_i*calc_k7_i*(calc_k4_i+calc_k5_i))+(calc_k4_i*calc_k6_i*(calc_k1_i+calc_k8_i)));	//211
    real calc_x3_i = ((calc_k1_i*calc_k3_i*(calc_k7_i+calc_k6_i))+(calc_k8_i*calc_k6_i*(calc_k2_i+calc_k3_i)));	//212
    real calc_x4_i = ((calc_k2_i*calc_k8_i*(calc_k4_i+calc_k5_i))+(calc_k3_i*calc_k5_i*(calc_k1_i+calc_k8_i)));	//213
    real calc_x1_ss = ((calc_k2_ss*calc_k4_ss*(calc_k7_ss+calc_k6_ss))+(calc_k5_ss*calc_k7_ss*(calc_k2_ss+calc_k3_ss)));	//247
    real calc_x2_ss = ((calc_k1_ss*calc_k7_ss*(calc_k4_ss+calc_k5_ss))+(calc_k4_ss*calc_k6_ss*(calc_k1_ss+calc_k8_ss)));	//248
    real calc_x3_ss = ((calc_k1_ss*calc_k3_ss*(calc_k7_ss+calc_k6_ss))+(calc_k8_ss*calc_k6_ss*(calc_k2_ss+calc_k3_ss)));	//249
    real calc_x4_ss = ((calc_k2_ss*calc_k8_ss*(calc_k4_ss+calc_k5_ss))+(calc_k3_ss*calc_k5_ss*(calc_k1_ss+calc_k8_ss)));	//250
    real calc_INaK = (calc_Pnak*((zna*calc_JnakNa)+(zk*calc_JnakK)));	//281
    real calc_E1_i = (calc_x1_i/(calc_x1_i+calc_x2_i+calc_x3_i+calc_x4_i));	//214
    real calc_E2_i = (calc_x2_i/(calc_x1_i+calc_x2_i+calc_x3_i+calc_x4_i));	//215
    real calc_E3_i = (calc_x3_i/(calc_x1_i+calc_x2_i+calc_x3_i+calc_x4_i));	//216
    real calc_E4_i = (calc_x4_i/(calc_x1_i+calc_x2_i+calc_x3_i+calc_x4_i));	//217
    real calc_E1_ss = (calc_x1_ss/(calc_x1_ss+calc_x2_ss+calc_x3_ss+calc_x4_ss));	//251
    real calc_E2_ss = (calc_x2_ss/(calc_x1_ss+calc_x2_ss+calc_x3_ss+calc_x4_ss));	//252
    real calc_E3_ss = (calc_x3_ss/(calc_x1_ss+calc_x2_ss+calc_x3_ss+calc_x4_ss));	//253
    real calc_E4_ss = (calc_x4_ss/(calc_x1_ss+calc_x2_ss+calc_x3_ss+calc_x4_ss));	//254
    real calc_JncxNa_i = (((3.000000000000000e+00f*((calc_E4_i*calc_k7_i)-(calc_E1_i*calc_k8_i)))+(calc_E3_i*calc_k4pp_i))-(calc_E2_i*calc_k3pp_i));	//219
    real calc_JncxCa_i = ((calc_E2_i*calc_k2_i)-(calc_E1_i*calc_k1_i));	//220
    real calc_JncxNa_ss = (((3.000000000000000e+00f*((calc_E4_ss*calc_k7_ss)-(calc_E1_ss*calc_k8_ss)))+(calc_E3_ss*calc_k4pp_ss))-(calc_E2_ss*calc_k3pp_ss));	//256
    real calc_JncxCa_ss = ((calc_E2_ss*calc_k2_ss)-(calc_E1_ss*calc_k1_ss));	//257
    real calc_INaCa_i = ((1.000000000000000e+00f-INaCa_fractionSS)*calc_Gncx*calc_allo_i*((zna*calc_JncxNa_i)+(zca*calc_JncxCa_i)));	//222
    real calc_INaCa_ss = (INaCa_fractionSS*calc_Gncx*calc_allo_ss*((zna*calc_JncxNa_ss)+(zca*calc_JncxCa_ss)));	//258

    rDY_[0] = (-(calc_INa+calc_INaL+calc_Ito+calc_ICaL+calc_ICaNa+calc_ICaK+calc_IKr+calc_IKs+calc_IK1+calc_INaCa_i+calc_INaCa_ss+calc_INaK+calc_INab+calc_IKb+calc_IpCa+calc_ICab+calc_IClCa+calc_IClb+calc_I_katp+stim_current));
    rDY_[1] = ((aCaMK*calc_CaMKb*(calc_CaMKb+CaMKt_old_))-(bCaMK*CaMKt_old_));
    rDY_[2] = ((((-(calc_INa+calc_INaL+(3.000000000000000e+00f*calc_INaCa_i)+calc_ICaNa_i+(3.000000000000000e+00f*calc_INaK)+calc_INab))*calc_Acap)/(F*calc_vmyo))+((calc_JdiffNa*calc_vss)/calc_vmyo));
    rDY_[3] = ((((-(calc_ICaNa_ss+(3.000000000000000e+00f*calc_INaCa_ss)))*calc_Acap)/(F*calc_vss))-calc_JdiffNa);
    rDY_[4] = ((((-(((calc_Ito+calc_IKr+calc_IKs+calc_IK1+calc_IKb+calc_I_katp+stim_current)-(2.000000000000000e+00f*calc_INaK))+calc_ICaK_i))*calc_Acap)/(F*calc_vmyo))+((calc_JdiffK*calc_vss)/calc_vmyo));
    rDY_[5] = ((((-calc_ICaK_ss)*calc_Acap)/(F*calc_vss))-calc_JdiffK);
    rDY_[6] = (calc_Bcai*(((((-((calc_ICaL_i+calc_IpCa+calc_ICab)-(2.000000000000000e+00f*calc_INaCa_i)))*calc_Acap)/(2.000000000000000e+00f*F*calc_vmyo))-((calc_Jup*calc_vnsr)/calc_vmyo))+((calc_Jdiff*calc_vss)/calc_vmyo)));
    rDY_[7] = (calc_Bcass*(((((-(calc_ICaL_ss-(2.000000000000000e+00f*calc_INaCa_ss)))*calc_Acap)/(2.000000000000000e+00f*F*calc_vss))+((calc_Jrel*calc_vjsr)/calc_vss))-calc_Jdiff));
    rDY_[8] = (calc_Jup-((calc_Jtr*calc_vjsr)/calc_vnsr));
    rDY_[9] = (calc_Bcajsr*(calc_Jtr-calc_Jrel));
    rDY_[10] = ((calc_mss-m_old_)/calc_tm);
    rDY_[11] = ((calc_hss-h_old_)/calc_th);
    rDY_[12] = ((calc_jss-j_old_)/calc_tj);
    rDY_[13] = ((calc_hssp-hp_old_)/calc_th);
    rDY_[14] = ((calc_jss-jp_old_)/calc_tjp);
    rDY_[15] = ((calc_mLss-mL_old_)/calc_tmL);
    rDY_[16] = ((calc_hLss-hL_old_)/thL);
    rDY_[17] = ((calc_hLssp-hLp_old_)/calc_thLp);
    rDY_[18] = ((calc_ass-a_old_)/calc_ta);
    rDY_[19] = ((calc_iss-iF_old_)/calc_tiF);
    rDY_[20] = ((calc_iss-iS_old_)/calc_tiS);
    rDY_[21] = ((calc_assp-ap_old_)/calc_ta);
    rDY_[22] = ((calc_iss-iFp_old_)/calc_tiFp);
    rDY_[23] = ((calc_iss-iSp_old_)/calc_tiSp);
    rDY_[24] = ((calc_dss-d_old_)/calc_td);
    rDY_[25] = ((calc_fss-ff_old_)/calc_tff);
    rDY_[26] = ((calc_fss-fs_old_)/calc_tfs);
    rDY_[27] = ((calc_fcass-fcaf_old_)/calc_tfcaf);
    rDY_[28] = ((calc_fcass-fcas_old_)/calc_tfcas);
    rDY_[29] = ((calc_jcass-jca_old_)/tjca);
    rDY_[30] = ((calc_fss-ffp_old_)/calc_tffp);
    rDY_[31] = ((calc_fcass-fcafp_old_)/calc_tfcafp);
    rDY_[32] = ((calc_anca_ss*k2n)-(nca_ss_old_*calc_km2n));
    rDY_[33] = ((calc_anca_i*k2n)-(nca_i_old_*calc_km2n));
    rDY_[34] = ((calc_beta*C2_old_)-(calc_alpha*C3_old_));
    rDY_[35] = (((calc_alpha*C3_old_)+(beta_1*C1_old_))-((calc_beta+alpha_1)*C2_old_));
    rDY_[36] = (((alpha_1*C2_old_)+(calc_beta_2*O_old_)+(calc_beta_ItoC2*I_old_))-((beta_1+calc_alpha_2+calc_alpha_C2ToI)*C1_old_));
    rDY_[37] = (((calc_alpha_2*C1_old_)+(calc_beta_i*I_old_))-((calc_beta_2+calc_alpha_i)*O_old_));
    rDY_[38] = (((calc_alpha_C2ToI*C1_old_)+(calc_alpha_i*O_old_))-((calc_beta_ItoC2+calc_beta_i)*I_old_));
    rDY_[39] = ((calc_xs1ss-xs1_old_)/calc_txs1);
    rDY_[40] = ((calc_xs2ss-xs2_old_)/calc_txs2);
    rDY_[41] = ((calc_Jrel_inf-Jrel_np_old_)/calc_tau_rel);
    rDY_[42] = ((calc_Jrel_infp-Jrel_p_old_)/calc_tau_relp);


}

inline __device__ void RHS_gpu_matlab(real *sv, real *rDY_, real stim_current, int threadID, real dt) {

    const int celltype = 0;
    real R = 8314.0f;
    real T = 310.0f;
    real F = 96485.0f;

    // cell geometry
    real L = 0.01;
    real rad = 0.0011;
    real vcell = 1000 * 3.14 * rad * rad * L;
    real Ageo = 2 * 3.14 * rad * rad + 2 * 3.14 * rad * L;
    real Acap = 2 * Ageo;
    real vmyo = 0.68 * vcell;
    real vnsr = 0.0552 * vcell;
    real vjsr = 0.0048 * vcell;
    real vss = 0.02 * vcell;

    real nao = 140.0f;
    real cao = 1.8f;
    real ko = 5.0f;


    real ICaL_fractionSS = 0.8f;
    real INaCa_fractionSS = 0.35f;


    real INa_Multiplier = 1.0f;
    real ICaL_Multiplier = 1.0f;
    real Ito_Multiplier = 1.0f;
    real INaL_Multiplier = 1.0f;
    real IKr_Multiplier = 1.0f;
    real IKs_Multiplier = 1.0f;
    real IK1_Multiplier = 1.0f;
    real IKb_Multiplier = 1.0f;
    real INaCa_Multiplier = 1.0f;
    real INaK_Multiplier = 1.0f;
    real INab_Multiplier = 1.0f;
    real ICab_Multiplier = 1.0f;
    real IpCa_Multiplier = 1.0f;
    real ICaCl_Multiplier = 1.0f;
    real IClb_Multiplier = 1.0f;
    real Jrel_Multiplier = 1.0f;
    real Jup_Multiplier = 1.0f;



    // give names to the state vector values
//    const real v       = *((real *)((char *)sv + pitch * 0) + threadID);
//    const real nai     = *((real *)((char *)sv + pitch * 1) + threadID);
//    const real nass    = *((real *)((char *)sv + pitch * 2) + threadID);
//    const real ki      = *((real *)((char *)sv + pitch * 3) + threadID);
//    const real kss     = *((real *)((char *)sv + pitch * 4) + threadID);
//    const real cai     = *((real *)((char *)sv + pitch * 5) + threadID);
//    const real cass    = *((real *)((char *)sv + pitch * 6) + threadID);
//    const real cansr   = *((real *)((char *)sv + pitch * 7) + threadID);
//    const real cajsr   = *((real *)((char *)sv + pitch * 8) + threadID);
//    const real m       = *((real *)((char *)sv + pitch * 9) + threadID);
//    const real hp      = *((real *)((char *)sv + pitch * 10) + threadID);
//    const real h       = *((real *)((char *)sv + pitch * 11) + threadID);
//    const real j       = *((real *)((char *)sv + pitch * 12) + threadID);
//    const real jp      = *((real *)((char *)sv + pitch * 13) + threadID);
//    const real mL      = *((real *)((char *)sv + pitch * 14) + threadID);
//    const real hL      = *((real *)((char *)sv + pitch * 15) + threadID);
//    const real hLp     = *((real *)((char *)sv + pitch * 16) + threadID);
//    const real a       = *((real *)((char *)sv + pitch * 17) + threadID);
//    const real iF      = *((real *)((char *)sv + pitch * 18) + threadID);
//    const real iS      = *((real *)((char *)sv + pitch * 19) + threadID);
//    const real ap      = *((real *)((char *)sv + pitch * 20) + threadID);
//    const real iFp     = *((real *)((char *)sv + pitch * 21) + threadID);
//    const real iSp     = *((real *)((char *)sv + pitch * 22) + threadID);
//    const real d       = *((real *)((char *)sv + pitch * 23) + threadID);
//    const real ff      = *((real *)((char *)sv + pitch * 24) + threadID);
//    const real fs      = *((real *)((char *)sv + pitch * 25) + threadID);
//    const real fcaf    = *((real *)((char *)sv + pitch * 26) + threadID);
//    const real fcas    = *((real *)((char *)sv + pitch * 27) + threadID);
//    const real jca     = *((real *)((char *)sv + pitch * 28) + threadID);
//    const real nca     = *((real *)((char *)sv + pitch * 29) + threadID);
//    const real nca_i   = *((real *)((char *)sv + pitch * 30) + threadID);
//    const real ffp     = *((real *)((char *)sv + pitch * 31) + threadID);
//    const real fcafp   = *((real *)((char *)sv + pitch * 32) + threadID);
//    const real xs1     = *((real *)((char *)sv + pitch * 33) + threadID);
//    const real xs2     = *((real *)((char *)sv + pitch * 34) + threadID);
//    const real Jrel_np = *((real *)((char *)sv + pitch * 35) + threadID);
//    const real CaMKt   = *((real *)((char *)sv + pitch * 36) + threadID);
//    const real ikr_c0  = *((real *)((char *)sv + pitch * 37) + threadID);
//    const real ikr_c1  = *((real *)((char *)sv + pitch * 38) + threadID);
//    const real ikr_c2  = *((real *)((char *)sv + pitch * 39) + threadID);
//    const real ikr_o   = *((real *)((char *)sv + pitch * 40) + threadID);
//    const real ikr_i   = *((real *)((char *)sv + pitch * 41) + threadID);
//    const real Jrel_p  = *((real *)((char *)sv + pitch * 42) + threadID);

    const real v       = sv[0];
    const real nai     = sv[1];
    const real nass    = sv[2];
    const real ki      = sv[3];
    const real kss     = sv[4];
    const real cai     = sv[5];
    const real cass    = sv[6];
    const real cansr   = sv[7];
    const real cajsr   = sv[8];
    const real m       = sv[9];
    const real hp      = sv[10];
    const real h       = sv[11];
    const real j       = sv[12];
    const real jp      = sv[13];
    const real mL      = sv[14];
    const real hL      = sv[15];
    const real hLp     = sv[16];
    const real a       = sv[17];
    const real iF      = sv[18];
    const real iS      = sv[19];
    const real ap      = sv[20];
    const real iFp     = sv[21];
    const real iSp     = sv[22];
    const real d       = sv[23];
    const real ff      = sv[24];
    const real fs      = sv[25];
    const real fcaf    = sv[26];
    const real fcas    = sv[27];
    const real jca     = sv[28];
    const real nca     = sv[29];
    const real nca_i   = sv[30];
    const real ffp     = sv[31];
    const real fcafp   = sv[32];
    const real xs1     = sv[33];
    const real xs2     = sv[34];
    const real Jrel_np = sv[35];
    const real CaMKt   = sv[36];
    const real ikr_c0  = sv[37];
    const real ikr_c1  = sv[38];
    const real ikr_c2  = sv[39];
    const real ikr_o   = sv[40];
    const real ikr_i   = sv[41];
    const real Jrel_p  = sv[42];

    real cli = 24;  // Intracellular Cl  [mM]
    real clo = 150; // Extracellular Cl  [mM]

    // CaMK constants
    real KmCaMK = 0.15;

    real aCaMK = 0.05;
    real bCaMK = 0.00068;
    real CaMKo = 0.05;
    real KmCaM = 0.0015;

    // update CaMK
    real CaMKb = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / cass);
    real CaMKa = CaMKb + CaMKt;
    real dCaMKt = aCaMK * CaMKb * (CaMKb + CaMKt) - bCaMK * CaMKt;

    // reversal potentials
    real ENa = (R * T / F) * log(nao / nai);
    real EK = (R * T / F) * log(ko / ki);
    real PKNa = 0.01833;
    real EKs = (R * T / F) * log((ko + PKNa * nao) / (ki + PKNa * nai));

    // convenient shorthand calculations
    real vffrt = v * F * F / (R * T);
    real vfrt = v * F / (R * T);

    real fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
    real fINaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
    real fItop = (1.0 / (1.0 + KmCaMK / CaMKa));
    real fICaLp = (1.0 / (1.0 + KmCaMK / CaMKa));

    // INa
    //[] = getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, INa_Multiplier);
    real INa, dm, dh, dhp, dj, djp;

    getINa_Grandi(&INa, &dm, &dh, &dhp, &dj, &djp, v, m, h, hp, j, jp, fINap, ENa, INa_Multiplier);

    //INaL
    real INaL, dmL, dhL, dhLp;
    //[INaL, dmL, dhL, dhLp] = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, celltype, INaL_Multiplier);

    getINaL_ORd2011(&INaL, &dmL, &dhL, &dhLp, v, mL, hL, hLp, fINaLp, ENa, celltype, INaL_Multiplier);

    //ITo
    // [Ito, da, diF, diS, dap, diFp, diSp] =
    // getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, Ito_Multiplier);

    real Ito, da, diF, diS, dap, diFp, diSp;
    getITo_ORd2011(&Ito, &da, &diF, &diS, &dap, &diFp, &diSp, v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, Ito_Multiplier);

    //ICaL
    //[ICaL_ss, ICaNa_ss, ICaK_ss, ICaL_i, ICaNa_i, ICaK_i, dd, dff, dfs, dfcaf, dfcas, djca, dnca, dnca_i,
    //         ... dffp, dfcafp, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] =
    //    getICaL_ORd2011_jt(v, d, ff, fs, fcaf, fcas, jca, nca, nca_i, ffp, fcafp, ... fICaLp, cai, cass, cao, nai, nass,
    //                       nao, ki, kss, ko, cli, clo, celltype, ICaL_fractionSS, ICaL_Multiplier);

    real ICaL_ss, ICaNa_ss, ICaK_ss, ICaL_i, ICaNa_i, ICaK_i, dd, dff, dfs, dfcaf, dfcas, djca, dnca, dnca_i, dffp, dfcafp, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo;

    getICaL_ORd2011_jt(&ICaL_ss, &ICaNa_ss, &ICaK_ss, &ICaL_i, &ICaNa_i, &ICaK_i, &dd, &dff, &dfs, &dfcaf, &dfcas, &djca, &dnca, &dnca_i,
                       &dffp, &dfcafp, &PhiCaL_ss, &PhiCaL_i, &gammaCaoMyo, &gammaCaiMyo, v, d, ff, fs, fcaf, fcas, jca, nca, nca_i, ffp, fcafp, fICaLp, cai, cass, cao, nai, nass,
                           nao, ki, kss, ko, cli, clo, celltype, ICaL_fractionSS, ICaL_Multiplier);

    real ICaL = ICaL_ss + ICaL_i;
    real ICaNa = ICaNa_ss + ICaNa_i;
    real ICaK = ICaK_ss + ICaK_i;

    //IKr[IKr, dt_ikr_c0, dt_ikr_c1, dt_ikr_c2, dt_ikr_o, dt_ikr_i] =
    //getIKr_ORd2011_MM(v, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i, ... ko, EK, celltype, IKr_Multiplier);
    real IKr, dikr_c0, dikr_c1, dikr_c2, dikr_o, dikr_i;
    getIKr_ORd2011_MM(&IKr, &dikr_c0, &dikr_c1, &dikr_c2, &dikr_o, &dikr_i, v, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i, ko, EK, celltype, IKr_Multiplier);

    //IKs
    // [IKs, dxs1, dxs2] = getIKs_ORd2011(v, xs1, xs2, cai, EKs, celltype, IKs_Multiplier);
    real IKs, dxs1, dxs2;
    getIKs_ORd2011(&IKs, &dxs1, &dxs2, v, xs1, xs2, cai, EKs, celltype, IKs_Multiplier);

    //IK1
    // IK1 = getIK1_CRLP(v, ko, EK, celltype, IK1_Multiplier);
    real IK1;
    getIK1_CRLP(&IK1, v, ko, EK, celltype, IK1_Multiplier);

    //INaCa[INaCa_i, INaCa_ss] =
    //getINaCa_ORd2011(v, F, R, T, nass, nai, nao, cass, cai, cao, celltype, INaCa_Multiplier, INaCa_fractionSS);

    real INaCa_i, INaCa_ss;
    getINaCa_ORd2011(&INaCa_i, &INaCa_ss, v, F, R, T, nass, nai, nao, cass, cai, cao, celltype, INaCa_Multiplier, INaCa_fractionSS);

    //INaK
    // INaK = getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, celltype, INaK_Multiplier);
    real INaK;
    getINaK_ORd2011(&INaK, v, F, R, T, nai, nao, ki, ko, celltype, INaK_Multiplier);

    //Minor background currents
    // calculate IKb
    real xkb = 1.0 / (1.0 + expf(-(v - 10.8968) / (23.9871)));
    real GKb = 0.0189 * IKb_Multiplier;
    if (celltype == 1)
        GKb = GKb * 0.6;

    real IKb = GKb * xkb * (v - EK);

    // calculate INab
    real PNab = 1.9239e-09 * INab_Multiplier;
    real INab = PNab * vffrt * (nai * expf(vfrt) - nao) / (expf(vfrt) - 1.0);

    // calculate ICab
    real PCab = 5.9194e-08 * ICab_Multiplier;
    //
    real ICab = PCab * 4.0 * vffrt * (gammaCaiMyo * cai * expf(2.0 * vfrt) - gammaCaoMyo * cao) / (expf(2.0 * vfrt) - 1.0);

    // calculate IpCa
    real GpCa = 5e-04 * IpCa_Multiplier;
    real IpCa = GpCa * cai / (0.0005 + cai);

    //Chloride
    // I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current

    real ecl = (R * T / F) * log(cli / clo); // [mV]

    real Fjunc = 1.0f;
    real Fsl = 1.0f - Fjunc; // fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace

    real GClCa = ICaCl_Multiplier * 0.2843; // [mS/uF]
    real GClB = IClb_Multiplier * 1.98e-3;  // [mS/uF] //
    real KdClCa = 0.1;                      // [mM]

    real I_ClCa_junc = Fjunc * GClCa / (1 + KdClCa / cass) * (v - ecl);
    real I_ClCa_sl = Fsl * GClCa / (1 + KdClCa / cai) * (v - ecl);

    real I_ClCa = I_ClCa_junc + I_ClCa_sl;
    real I_Clbk = GClB * (v - ecl);

    //Calcium handling
    // calculate ryanodione receptor calcium induced calcium release from the jsr
    real fJrelp = (1.0 / (1.0 + KmCaMK / CaMKa));

    //Jrel[Jrel, dJrel_np, dJrel_p] =
    //getJrel_ORd2011(Jrel_np, Jrel_p, ICaL_ss, cass, cajsr, fJrelp, celltype, Jrel_Multiplier);
    real Jrel, dJrel_np, dJrel_p;
    getJrel_ORd2011(&Jrel, &dJrel_np, &dJrel_p, Jrel_np, Jrel_p, ICaL_ss, cass, cajsr, fJrelp, celltype, Jrel_Multiplier);

    real fJupp = (1.0 / (1.0 + KmCaMK / CaMKa));
    //[ Jup, Jleak ] = getJup_ORd2011(cai, cansr, fJupp, celltype, Jup_Multiplier);
    real Jup, Jleak;
    getJup_ORd2011(&Jup, &Jleak, cai, cansr, fJupp, celltype, Jup_Multiplier);

    // calculate tranlocation flux
    real Jtr = (cansr - cajsr) / 60.0f;
    // update the membrane voltage

    real dv = -(INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa_i + INaCa_ss + INaK + INab + IKb + IpCa + ICab + +I_ClCa + I_Clbk + stim_current);

    real JdiffNa = (nass - nai) / 2.0f;
    real JdiffK = (kss - ki) / 2.0f;
    real Jdiff = (cass - cai) / 0.2f;

    // calcium buffer constants
    real cmdnmax = 0.05;
    if (celltype == 1)
        cmdnmax = cmdnmax * 1.3;

    real kmcmdn = 0.00238;
    real trpnmax = 0.07;
    real kmtrpn = 0.0005;
    real BSRmax = 0.047;
    real KmBSR = 0.00087;
    real BSLmax = 1.124;
    real KmBSL = 0.0087;
    real csqnmax = 10.0;
    real kmcsqn = 0.8;

    // update intracellular concentrations, using buffers for cai, cass, cajsr
    real dnai = -(ICaNa_i + INa + INaL + 3.0 * INaCa_i + 3.0 * INaK + INab) * Acap / (F * vmyo) + JdiffNa * vss / vmyo;
    real dnass = -(ICaNa_ss + 3.0 * INaCa_ss) * Acap / (F * vss) - JdiffNa;

    real dki = -(ICaK_i + Ito + IKr + IKs + IK1 + IKb + stim_current - 2.0 * INaK) * Acap / (F * vmyo) + JdiffK * vss / vmyo;
    real dkss = -(ICaK_ss)*Acap / (F * vss) - JdiffK;

    real Bcai = 1.0 / (1.0 + cmdnmax * kmcmdn / powf((kmcmdn + cai), 2.0f) + trpnmax * kmtrpn / powf((kmtrpn + cai), 2.0f));
    real dcai = Bcai * (-(ICaL_i + IpCa + ICab - 2.0 * INaCa_i) * Acap / (2.0 * F * vmyo) - Jup * vnsr / vmyo + Jdiff * vss / vmyo);

    real Bcass = 1.0 / (1.0 + BSRmax * KmBSR / powf((KmBSR + cass), 2.0f) + BSLmax * KmBSL / powf((KmBSL + cass), 2.0f));
    real dcass = Bcass * (-(ICaL_ss - 2.0 * INaCa_ss) * Acap / (2.0 * F * vss) + Jrel * vjsr / vss - Jdiff);

    real dcansr = Jup - Jtr * vjsr / vnsr;

    real Bcajsr = 1.0 / (1.0 + csqnmax * kmcsqn / powf((kmcsqn + cajsr), 2.0f));
    real dcajsr = Bcajsr * (Jtr - Jrel);

    rDY_[0]  = dv      ;
    rDY_[1]  = dnai    ;
    rDY_[2]  = dnass   ;
    rDY_[3]  = dki     ;
    rDY_[4]  = dkss    ;
    rDY_[5]  = dcai    ;
    rDY_[6]  = dcass   ;
    rDY_[7]  = dcansr  ;
    rDY_[8]  = dcajsr  ;
    rDY_[9]  = dm      ;
    rDY_[10] = dhp     ;
    rDY_[11] = dh      ;
    rDY_[12] = dj      ;
    rDY_[13] = djp     ;
    rDY_[14] = dmL     ;
    rDY_[15] = dhL     ;
    rDY_[16] = dhLp    ;
    rDY_[17] = da      ;
    rDY_[18] = diF     ;
    rDY_[19] = diS     ;
    rDY_[20] = dap     ;
    rDY_[21] = diFp    ;
    rDY_[22] = diSp    ;
    rDY_[23] = dd      ;
    rDY_[24] = dff     ;
    rDY_[25] = dfs     ;
    rDY_[26] = dfcaf   ;
    rDY_[27] = dfcas   ;
    rDY_[28] = djca    ;
    rDY_[29] = dnca    ;
    rDY_[30] = dnca_i  ;
    rDY_[31] = dffp    ;
    rDY_[32] = dfcafp  ;
    rDY_[33] = dxs1    ;
    rDY_[34] = dxs2    ;
    rDY_[35] = dJrel_np;
    rDY_[36] = dCaMKt  ;
    rDY_[37] = dikr_c0 ;
    rDY_[38] = dikr_c1 ;
    rDY_[39] = dikr_c2 ;
    rDY_[40] = dikr_o  ;
    rDY_[41] = dikr_i  ;
    rDY_[42] = dJrel_p ;

}
    //INa formulations
__inline__ __device__ void getINa_Grandi(real *INa, real *dm, real *dh, real *dhp, real *dj, real *djp, real v, real m, real h, real hp, real j, real jp, real fINap, real ENa, real INa_Multiplier) {

    // The Grandi implementation updated with INa phosphorylation.
    //m gate

    real mss = 1.0f / powf((1 + expf(-(56.86 + v) / 9.03)), 2.0f);
    real taum = 0.1292 * expf(-powf(((v + 45.79f) / 15.54f), 2.0f)) + 0.06487 * expf(-powf(((v - 4.823) / 51.12), 2.0f));
    *dm = (mss - m) / taum;

    //h gate
    real ah = (v >= -40) ? 0 : (0.057 * expf(-(v + 80) / 6.8));

    real bh = (v >= -40) ? (0.77 / (0.13 * (1 + expf(-(v + 10.66) / 11.1)))) : ((2.7 * expf(0.079 * v) + 3.1e+05 * expf(0.3485 * v)));

    real tauh = 1.0f / (ah + bh);
    real hss = 1.0f / powf( (1.0f + expf((v + 71.55) / 7.43)), 2.0f);
    *dh = (hss - h) / tauh;
    //j gate
    real aj = (v >= -40) ? 0 : (((-2.5428e+04 * expf(0.2444 * v) - 6.948e-06 * expf(-0.04391 * v)) * (v + 37.78))/(1 + expf(0.311 * (v + 79.23))));
    real bj = (v >= -40) ? ((0.6 * expf(0.057 * v)) / (1 + expf(-0.1 * (v + 32)))) :((0.02424 * expf(-0.01052 * v)) / (1 + expf(-0.1378 * (v + 40.14))));
    real tauj = 1.0f / (aj + bj);
    real jss = 1.0f / powf((1 + expf((v + 71.55) / 7.43)), 2.0f);
    *dj = (jss - j) / tauj;

    //h phosphorylated
    real hssp = 1.0f / ( powf((1.0f + expf((v + 71.55f + 6.0f) / 7.43)), 2.0f) );
    *dhp = (hssp - hp) / tauh;
    //j phosphorylated
    real taujp = 1.46f * tauj;
    *djp = (jss - jp) / taujp;

    real GNa = 11.7802f;
    *INa = INa_Multiplier * GNa * (v - ENa) * powf(m, 3.0f) * ((1.0 - fINap) * h * j + fINap * hp * jp);
}

//INaL
__inline__ __device__ void getINaL_ORd2011(real *INaL, real *dmL, real *dhL, real *dhLp, real v, real mL, real hL, real hLp, real fINaLp, real ENa, int celltype, real INaL_Multiplier) {
    // calculate INaL
    real mLss = 1.0 / (1.0 + expf((-(v + 42.85)) / 5.264));
    real tm = 0.1292 * expf(-powf(((v + 45.79) / 15.54), 2)) + 0.06487 * expf(-powf(((v - 4.823) / 51.12), 2.0f));
    real tmL = tm;
    *dmL = (mLss - mL) / tmL;
    real hLss = 1.0 / (1.0 + expf((v + 87.61) / 7.488));
    real thL = 200.0;
    *dhL = (hLss - hL) / thL;
    real hLssp = 1.0 / (1.0 + expf((v + 93.81) / 7.488));
    real thLp = 3.0 * thL;
    *dhLp = (hLssp - hLp) / thLp;

    real GNaL = 0.0279 * INaL_Multiplier;

    if (celltype == 1) GNaL = GNaL * 0.6;

    *INaL = GNaL * (v - ENa) * mL * ((1.0 - fINaLp) * hL + fINaLp * hLp);
}

//ITo
__inline__ __device__ void getITo_ORd2011(real *Ito, real *da, real *diF, real *diS, real *dap, real *diFp, real *diSp, real v, real a,
                    real iF, real iS, real ap, real iFp, real iSp, real fItop, real EK, int celltype, real Ito_Multiplier) {

    // calculate Ito
    real ass = 1.0 / (1.0 + expf((-(v - 14.34)) / 14.82));
    real ta =
        1.0515 / (1.0 / (1.2089 * (1.0 + expf(-(v - 18.4099) / 29.3814))) + 3.5 / (1.0 + expf((v + 100.0) / 29.3814)));
    *da = (ass - a) / ta;
    real iss = 1.0 / (1.0 + expf((v + 43.94) / 5.711));

    real delta_epi;

    if(celltype == 1)
        delta_epi = 1.0 - (0.95 / (1.0 + expf((v + 70.0) / 5.0)));
    else
        delta_epi = 1.0;

    real tiF = 4.562 + 1 / (0.3933 * expf((-(v + 100.0)) / 100.0) + 0.08004 * expf((v + 50.0) / 16.59));
    real tiS = 23.62 + 1 / (0.001416 * expf((-(v + 96.52)) / 59.05) + 1.780e-8 * expf((v + 114.1) / 8.079));

    tiF = tiF * delta_epi;
    tiS = tiS * delta_epi;

    real AiF = 1.0 / (1.0 + expf((v - 213.6) / 151.2));
    real AiS = 1.0 - AiF;
    *diF = (iss - iF) / tiF;
    *diS = (iss - iS) / tiS;
    real i = AiF * iF + AiS * iS;
    real assp = 1.0 / (1.0 + expf((-(v - 24.34)) / 14.82));
    *dap = (assp - ap) / ta;
    real dti_develop = 1.354 + 1.0e-4 / (expf((v - 167.4) / 15.89) + expf(-(v - 12.23) / 0.2154));
    real dti_recover = 1.0 - 0.5 / (1.0 + expf((v + 70.0) / 20.0));
    real tiFp = dti_develop * dti_recover * tiF;
    real tiSp = dti_develop * dti_recover * tiS;

    *diFp = (iss - iFp) / tiFp;
    *diSp = (iss - iSp) / tiSp;
    real ip = AiF * iFp + AiS * iSp;
    real Gto = 0.16 * Ito_Multiplier;

    if(celltype == 1)
        Gto = Gto * 2.0;
    else if(celltype == 2)
        Gto = Gto * 2.0;

    *Ito = Gto * (v - EK) * ((1.0 - fItop) * a * i + fItop * ap * ip);
}


    // a variant updated by jakub, using a changed activation curve
                                                             // it computes both ICaL in subspace and myoplasm (_i)
__inline__ __device__ void getICaL_ORd2011_jt(real *ICaL_ss, real *ICaNa_ss, real *ICaK_ss, real *ICaL_i, real *ICaNa_i, real *ICaK_i,
                        real *dd, real *dff, real *dfs, real *dfcaf,real * dfcas,real * djca, real *dnca, real * dnca_i,
                        real *dffp, real *dfcafp, real *PhiCaL_ss, real *PhiCaL_i, real *gammaCaoMyo, real *gammaCaiMyo,
                        real v, real d, real ff, real fs, real fcaf, real fcas, real jca, real nca, real nca_i, real ffp,
                        real fcafp, real fICaLp, real cai, real cass, real cao, real nai, real nass, real nao, real ki,
                        real kss, real ko, real cli, real clo, real celltype, real ICaL_fractionSS, real ICaL_PCaMultiplier) {

    // physical constants
    real R = 8314.0;
    real T = 310.0;
    real F = 96485.0;
    real vffrt = v * F * F / (R * T);
    real vfrt = v * F / (R * T);

    // calculate ICaL, ICaNa, ICaK

    real dss = 1.0763 * expf(-1.0070 * expf(-0.0829 * (v))); // magyar
    if(v > 31.4978)                                        // activation cannot be greater than 1
        dss = 1.0f;

    real td = 0.6 + 1.0 / (expf(-0.05 * (v + 6.0)) + expf(0.09 * (v + 14.0)));

    *dd = (dss - d) / td;
    real fss = 1.0 / (1.0 + expf((v + 19.58) / 3.696));
    real tff = 7.0 + 1.0 / (0.0045 * expf(-(v + 20.0) / 10.0) + 0.0045 * expf((v + 20.0) / 10.0));
    real tfs = 1000.0 + 1.0 / (0.000035 * expf(-(v + 5.0) / 4.0) + 0.000035 * expf((v + 5.0) / 6.0));
    real Aff = 0.6;
    real Afs = 1.0 - Aff;
    *dff = (fss - ff) / tff;
    *dfs = (fss - fs) / tfs;
    real f = Aff * ff + Afs * fs;
    real fcass = fss;
    real tfcaf = 7.0 + 1.0 / (0.04 * expf(-(v - 4.0) / 7.0) + 0.04 * expf((v - 4.0) / 7.0));
    real tfcas = 100.0 + 1.0 / (0.00012 * expf(-v / 3.0) + 0.00012 * expf(v / 7.0));

    real Afcaf = 0.3 + 0.6 / (1.0 + expf((v - 10.0) / 10.0));

    real Afcas = 1.0 - Afcaf;
    *dfcaf = (fcass - fcaf) / tfcaf;
    *dfcas = (fcass - fcas) / tfcas;
    real fca = Afcaf * fcaf + Afcas * fcas;

    real tjca = 75;
    real jcass = 1.0 / (1.0 + expf((v + 18.08) / (2.7916)));
    *djca = (jcass - jca) / tjca;
    real tffp = 2.5 * tff;
    *dffp = (fss - ffp) / tffp;
    real fp = Aff * ffp + Afs * fs;
    real tfcafp = 2.5 * tfcaf;
    *dfcafp = (fcass - fcafp) / tfcafp;
    real fcap = Afcaf * fcafp + Afcas * fcas;

    // SS nca
    real Kmn = 0.002;
    real k2n = 500.0;
    real km2n = jca * 1;
    real anca = 1.0 / (k2n / km2n + powf((1.0 + Kmn / cass), 4.0f));
    *dnca = anca * k2n - nca * km2n;

    // myoplasmic nca
    real anca_i = 1.0 / (k2n / km2n + powf((1.0 + Kmn / cai), 4.0f));
    *dnca_i = anca_i * k2n - nca_i * km2n;

    // SS driving force
    clo = 150;
    cli = 24;
    real Io =
        0.5 * (nao + ko + clo + 4 * cao) / 1000; // ionic strength outside. /1000 is for things being in micromolar
    real Ii =
        0.5 * (nass + kss + cli + 4 * cass) / 1000; // ionic strength outside. /1000 is for things being in micromolar
                                                    // The ionic strength is too high for basic DebHuc. We'll use Davies
    real dielConstant = 74;                         // water at 37.
    real temp = 310;                                // body temp in kelvins.
    real constA = 1.82e+06f * powf((dielConstant * temp), -1.5f);

    real gamma_cai = expf(-constA * 4 * (sqrtf(Ii) / (1.0f + sqrtf(Ii)) - 0.3f * Ii));
    real gamma_cao = expf(-constA * 4 * (sqrtf(Io) / (1.0f + sqrtf(Io)) - 0.3f * Io));
    real gamma_nai = expf(-constA * (sqrtf(Ii) / (1.0f + sqrtf(Ii)) - 0.3f * Ii));
    real gamma_nao = expf(-constA * (sqrtf(Io) / (1.0f + sqrtf(Io)) - 0.3f * Io));
    real gamma_ki  = expf(-constA *  (sqrtf(Ii) / (1.0f + sqrtf(Ii)) - 0.3f * Ii));
    real gamma_kao = expf(-constA * (sqrtf(Io) / (1.0f + sqrtf(Io)) - 0.3f * Io));

    *PhiCaL_ss = 4.0 * vffrt * (gamma_cai * cass * expf(2.0 * vfrt) - gamma_cao * cao) / (expf(2.0 * vfrt) - 1.0);
    real PhiCaNa_ss = 1.0 * vffrt * (gamma_nai * nass * expf(1.0 * vfrt) - gamma_nao * nao) / (expf(1.0 * vfrt) - 1.0);
    real PhiCaK_ss = 1.0 * vffrt * (gamma_ki * kss * expf(1.0 * vfrt) - gamma_kao * ko) / (expf(1.0 * vfrt) - 1.0);

    // Myo driving force
    Io = 0.5 * (nao + ko + clo + 4 * cao) / 1000; // ionic strength outside. /1000 is for things being in micromolar
    Ii = 0.5 * (nai + ki + cli + 4 * cai) / 1000; // ionic strength outside. /1000 is for things being in micromolar
                                                  // The ionic strength is too high for basic DebHuc. We'll use Davies
    dielConstant = 74.0f;                            // water at 37
    temp = 310.0f;                                   // body temp in kelvins.
    constA = 1.82e+06f * powf((dielConstant * temp), -1.5f);

    gamma_cai = expf(-constA * 4 * (sqrtf(Ii) / (1 + sqrtf(Ii)) - 0.3 * Ii));
    gamma_cao = expf(-constA * 4 * (sqrtf(Io) / (1 + sqrtf(Io)) - 0.3 * Io));
    gamma_nai = expf(-constA * (sqrtf(Ii) / (1 + sqrtf(Ii)) - 0.3 * Ii));
    gamma_nao = expf(-constA * (sqrtf(Io) / (1 + sqrtf(Io)) - 0.3 * Io));
    gamma_ki = expf(-constA  * (sqrtf(Ii) / (1 + sqrtf(Ii)) - 0.3 * Ii));
    gamma_kao = expf(-constA * (sqrtf(Io) / (1 + sqrtf(Io)) - 0.3 * Io));

    *gammaCaoMyo = gamma_cao;
    *gammaCaiMyo = gamma_cai;

    *PhiCaL_i = 4.0 * vffrt * (gamma_cai * cai * expf(2.0 * vfrt) - gamma_cao * cao) / (expf(2.0 * vfrt) - 1.0);
    real PhiCaNa_i = 1.0 * vffrt * (gamma_nai * nai * expf(1.0 * vfrt) - gamma_nao * nao) / (expf(1.0 * vfrt) - 1.0);
    real PhiCaK_i = 1.0 * vffrt * (gamma_ki * ki * expf(1.0 * vfrt) - gamma_kao * ko) / (expf(1.0 * vfrt) - 1.0);
    // The rest
    real PCa = 8.3757e-05 * ICaL_PCaMultiplier;

    if(celltype == 1)
        PCa = PCa * 1.2f;
    else if(celltype == 2)
        PCa = PCa * 2.0f;

    real PCap = 1.1 * PCa;
    real PCaNa = 0.00125 * PCa;
    real PCaK = 3.574e-4 * PCa;
    real PCaNap = 0.00125 * PCap;
    real PCaKp = 3.574e-4 * PCap;

    *ICaL_ss = (1.0f - fICaLp) * PCa * (*PhiCaL_ss) * d * (f * (1.0 - nca) + jca * fca * nca) +
               fICaLp * PCap * (*PhiCaL_ss) * d * (fp * (1.0 - nca) + jca * fcap * nca);

    *ICaNa_ss = (1.0 - fICaLp) * PCaNa * PhiCaNa_ss * d * (f * (1.0 - nca) + jca * fca * nca) +
                fICaLp * PCaNap * PhiCaNa_ss * d * (fp * (1.0 - nca) + jca * fcap * nca);

    *ICaK_ss = (1.0 - fICaLp) * PCaK * PhiCaK_ss * d * (f * (1.0 - nca) + jca * fca * nca) +
               fICaLp * PCaKp * PhiCaK_ss * d * (fp * (1.0 - nca) + jca * fcap * nca);

    *ICaL_i = (1.0 - fICaLp) * PCa * (*PhiCaL_i) * d * (f * (1.0 - nca_i) + jca * fca * nca_i) +
             fICaLp * PCap * (*PhiCaL_i) * d * (fp * (1.0 - nca_i) + jca * fcap * nca_i);

    *ICaNa_i = (1.0 - fICaLp) * PCaNa * PhiCaNa_i * d * (f * (1.0 - nca_i) + jca * fca * nca_i) +
               fICaLp * PCaNap * PhiCaNa_i * d * (fp * (1.0 - nca_i) + jca * fcap * nca_i);

    *ICaK_i = (1.0 - fICaLp) * PCaK * PhiCaK_i * d * (f * (1.0 - nca_i) + jca * fca * nca_i) +
              fICaLp * PCaKp * PhiCaK_i * d * (fp * (1.0 - nca_i) + jca * fcap * nca_i);

    // And we weight ICaL (in ss) and ICaL_i
    *ICaL_i = *ICaL_i * (1 - ICaL_fractionSS);
    *ICaNa_i = *ICaNa_i * (1 - ICaL_fractionSS);
    *ICaK_i = *ICaK_i * (1 - ICaL_fractionSS);
    *ICaL_ss = *ICaL_ss * ICaL_fractionSS;
    *ICaNa_ss = *ICaNa_ss * ICaL_fractionSS;
    *ICaK_ss = *ICaK_ss * ICaL_fractionSS;
}

    // Variant based on Lu-Vandenberg
__inline__ __device__ void getIKr_ORd2011_MM(real *IKr, real *dc0, real *dc1, real *dc2, real *do_, real *di, real V, real c0,
                                    real c1, real c2, real  o, real i, real ko, real EK, int celltype, real IKr_Multiplier) {

    // physical constants
    real R = 8314.0;
    real T = 310.0;
    real F = 96485.0;

    real vfrt = V * F / (R * T);

    // transition rates
    // from c0 to c1 in l-v model,
    real alpha = 0.1161 * expf(0.2990 * vfrt);
    // from c1 to c0 in l-v/
    real beta = 0.2442 * expf(-1.604 * vfrt);

    // from c1 to c2 in l-v/
    real alpha1 = 1.25 * 0.1235;
    // from c2 to c1 in l-v/
    real beta1 = 0.1911;

    // from c2 to o/           c1 to o
    real alpha2 = 0.0578 * expf(0.9710 * vfrt); //
    // from o to c2/
    real beta2 = 0.349e-3 * expf(-1.062 * vfrt); //

    // from o to i
    real alphai = 0.2533 * expf(0.5953 * vfrt); //
    // from i to o
    real betai = 1.25 * 0.0522 * expf(-0.8209 * vfrt); //

    // from c2 to i (from c1 in orig)
    real alphac2ToI = 0.52e-4 * expf(1.525 * vfrt); //
    // from i to c2
    // betaItoC2 = 0.85e-8 * expf(-1.842 * vfrt); //
    real betaItoC2 = (beta2 * betai * alphac2ToI) / (alpha2 * alphai); //
    // transitions themselves
    // for reason of backward compatibility of naming of an older version of a MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

    *dc0 = c1 * beta - c0 * alpha;                         // delta for c0
    *dc1 = c0 * alpha + c2 * beta1 - c1 * (beta + alpha1); // c1
    *dc2 = c1 * alpha1 + o * beta2 + i * betaItoC2 -
           c2 * (beta1 + alpha2 + alphac2ToI); // subtraction is into c2, to o, to i. // c2
    *do_ = c2 * alpha2 + i * betai - o * (beta2 + alphai);
    *di = c2 * alphac2ToI + o * alphai - i * (betaItoC2 + betai);

    real GKr = 0.0321 * sqrtf(ko / 5) * IKr_Multiplier; // 1st element compensates for change to ko (sqrtf(5/5.4)* 0.0362)

    if(celltype == 1)
        GKr = GKr * 1.3;
    else if(celltype == 2)
        GKr = GKr * 0.8;

    *IKr = GKr * o * (V - EK);
}

__inline__ __device__ void getIKs_ORd2011(real *IKs, real *dxs1, real *dxs2, real v, real xs1, real xs2, real cai, real EKs, int celltype, real IKs_Multiplier) {
    // calculate IKs
    real xs1ss = 1.0 / (1.0 + expf((-(v + 11.60)) / 8.932));
    real txs1 = 817.3 + 1.0 / (2.326e-4 * expf((v + 48.28) / 17.80) + 0.001292 * expf((-(v + 210.0)) / 230.0));
    *dxs1 = (xs1ss - xs1) / txs1;
    real xs2ss = xs1ss;
    real txs2 = 1.0 / (0.01 * expf((v - 50.0) / 20.0) + 0.0193 * expf((-(v + 66.54)) / 31.0));
    *dxs2 = (xs2ss - xs2) / txs2;
    real KsCa = 1.0 + 0.6 / (1.0 + powf((3.8e-5 / cai), 1.4));
    real GKs = 0.0011 * IKs_Multiplier;
    if(celltype == 1)
        GKs = GKs * 1.4;
    *IKs = GKs * KsCa * xs1 * xs2 * (v - EKs);
}

__inline__ __device__ void getIK1_CRLP(real *IK1, real v, real ko, real  EK, int celltype, real IK1_Multiplier) {
                         // IK1
    real aK1 = 4.094/(1+expf(0.1217*(v-EK-49.934)));
    real bK1 = (15.72*expf(0.0674*(v-EK-3.257))+expf(0.0618*(v-EK-594.31)))/(1+expf(-0.1629*(v-EK+14.207)));
    real K1ss = aK1/(aK1+bK1);

    real GK1 = IK1_Multiplier  * 0.6992; //0.7266; //* sqrtf(5/5.4))
    if (celltype==1)
        GK1=GK1*1.2;
    else if (celltype==2)
        GK1=GK1*1.3;

    *IK1=GK1*sqrtf(ko/5)*K1ss*(v-EK);
}

__inline__ __device__ void getINaCa_ORd2011(real *INaCa_i, real *INaCa_ss, real v, real F, real R, real T, real nass,
                                            real nai, real nao, real cass, real cai, real cao, int celltype,
                                            real INaCa_Multiplier, real INaCa_fractionSS) {
    real zca = 2.0f;
    real kna1 = 15.0f;
    real kna2 = 5.0f;
    real kna3 = 88.12f;
    real kasymm = 12.5f;
    real wna = 6.0e4f;
    real wca = 6.0e4f;
    real wnaca = 5.0e3f;
    real kcaon = 1.5e6f;
    real kcaoff = 5.0e3f;
    real qna = 0.5224f;
    real qca = 0.1670f;
    real hca = expf((qca * v * F) / (R * T));
    real hna = expf((qna * v * F) / (R * T));
    real h1 = 1 + nai / kna3 * (1 + hna);
    real h2 = (nai * hna) / (kna3 * h1);
    real h3 = 1.0f / h1;
    real h4 = 1.0f + nai / kna1 * (1 + nai / kna2);
    real h5 = nai * nai / (h4 * kna1 * kna2);
    real h6 = 1.0f / h4;
    real h7 = 1.0f + nao / kna3 * (1.0f + 1.0f / hna);
    real h8 = nao / (kna3 * hna * h7);
    real h9 = 1.0f / h7;
    real h10 = kasymm + 1.0f + nao / kna1 * (1.0f + nao / kna2);
    real h11 = nao * nao / (h10 * kna1 * kna2);
    real h12 = 1.0 / h10;
    real k1 = h12 * cao * kcaon;
    real k2 = kcaoff;
    real k3p = h9 * wca;
    real k3pp = h8 * wnaca;
    real k3 = k3p + k3pp;
    real k4p = h3 * wca / hca;
    real k4pp = h2 * wnaca;
    real k4 = k4p + k4pp;
    real k5 = kcaoff;
    real k6 = h6 * cai * kcaon;
    real k7 = h5 * h2 * wna;
    real k8 = h8 * h11 * wna;
    real x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
    real x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
    real x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
    real x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
    real E1 = x1 / (x1 + x2 + x3 + x4);
    real E2 = x2 / (x1 + x2 + x3 + x4);
    real E3 = x3 / (x1 + x2 + x3 + x4);
    real E4 = x4 / (x1 + x2 + x3 + x4);
    real KmCaAct = 150.0e-6;
    real allo = 1.0 / (1.0 + powf((KmCaAct / cai), 2.0));
    real zna = 1.0;
    real JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
    real JncxCa = E2 * k2 - E1 * k1;
    real Gncx = 0.0034 * INaCa_Multiplier;
    if(celltype == 1)
        Gncx = Gncx * 1.1;
    else if(celltype == 2)
        Gncx = Gncx * 1.4;

    *INaCa_i = (1 - INaCa_fractionSS) * Gncx * allo * (zna * JncxNa + zca * JncxCa);

    // calculate INaCa_ss
    h1 = 1 + nass / kna3 * (1 + hna);
    h2 = (nass * hna) / (kna3 * h1);
    h3 = 1.0 / h1;
    h4 = 1.0 + nass / kna1 * (1 + nass / kna2);
    h5 = nass * nass / (h4 * kna1 * kna2);
    h6 = 1.0 / h4;
    h7 = 1.0 + nao / kna3 * (1.0 + 1.0 / hna);
    h8 = nao / (kna3 * hna * h7);
    h9 = 1.0 / h7;
    h10 = kasymm + 1.0 + nao / kna1 * (1 + nao / kna2);
    h11 = nao * nao / (h10 * kna1 * kna2);
    h12 = 1.0 / h10;
    k1 = h12 * cao * kcaon;
    k2 = kcaoff;
    k3p = h9 * wca;
    k3pp = h8 * wnaca;
    k3 = k3p + k3pp;
    k4p = h3 * wca / hca;
    k4pp = h2 * wnaca;
    k4 = k4p + k4pp;
    k5 = kcaoff;
    k6 = h6 * cass * kcaon;
    k7 = h5 * h2 * wna;
    k8 = h8 * h11 * wna;
    x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
    x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
    x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
    x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
    E1 = x1 / (x1 + x2 + x3 + x4);
    E2 = x2 / (x1 + x2 + x3 + x4);
    E3 = x3 / (x1 + x2 + x3 + x4);
    E4 = x4 / (x1 + x2 + x3 + x4);
    KmCaAct = 150.0e-6;
    allo = 1.0 / (1.0 + powf((KmCaAct / cass), 2.0f));
    JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
    JncxCa = E2 * k2 - E1 * k1;
    *INaCa_ss = INaCa_fractionSS * Gncx * allo * (zna * JncxNa + zca * JncxCa);
}


__inline__ __device__ void getINaK_ORd2011(real *INaK, real v, real F, real R, real T, real nai, real nao, real ki, real ko, int celltype, real  INaK_Multiplier) {

    // calculate INaK
    real zna = 1.0;
    real k1p = 949.5;
    real k1m = 182.4;
    real k2p = 687.2;
    real k2m = 39.4;
    real k3p = 1899.0;
    real k3m = 79300.0;
    real k4p = 639.0;
    real k4m = 40.0;
    real Knai0 = 9.073;
    real Knao0 = 27.78;
    real delta = -0.1550;
    real Knai = Knai0 * expf((delta * v * F) / (3.0 * R * T));
    real Knao = Knao0 * expf(((1.0 - delta) * v * F) / (3.0 * R * T));
    real Kki = 0.5;
    real Kko = 0.3582;
    real MgADP = 0.05;
    real MgATP = 9.8;
    real Kmgatp = 1.698e-7;
    real H = 1.0e-7;
    real eP = 4.2;
    real Khp = 1.698e-7;
    real Knap = 224.0;
    real Kxkur = 292.0;
    real P = eP / (1.0 + H / Khp + nai / Knap + ki / Kxkur);
    real a1 = (k1p *  powf((nai / Knai), 3.0f) ) / ( powf((1.0 + nai / Knai), 3.0f) + powf((1.0 + ki / Kki), 2.0f) - 1.0);
    real b1 = k1m * MgADP;
    real a2 = k2p;
    real b2 = (k2m * powf((nao / Knao), 3.0)) / ( powf((1.0 + nao / Knao), 3.0f) + powf((1.0 + ko / Kko), 2.0f) - 1.0);
    real a3 = (k3p * powf((ko / Kko), 2.0)) / ( powf((1.0 + nao / Knao), 3.0f) + powf((1.0 + ko / Kko), 2.0) - 1.0);
    real b3 = (k3m * P * H) / (1.0 + MgATP / Kmgatp);
    real a4 = (k4p * MgATP / Kmgatp) / (1.0 + MgATP / Kmgatp);
    real b4 = (k4m * powf((ki / Kki), 2.0f)) / ( powf((1.0 + nai / Knai), 3.0f) + powf((1.0f + ki / Kki), 2.0f) - 1.0);
    real x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2;
    real x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4;
    real x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1;
    real x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1;
    real E1 = x1 / (x1 + x2 + x3 + x4);
    real E2 = x2 / (x1 + x2 + x3 + x4);
    real E3 = x3 / (x1 + x2 + x3 + x4);
    real E4 = x4 / (x1 + x2 + x3 + x4);
    real zk = 1.0;
    real JnakNa = 3.0 * (E1 * a3 - E2 * b3);
    real JnakK = 2.0 * (E4 * b1 - E3 * a1);
    real Pnak = 15.4509 * INaK_Multiplier;

    if(celltype == 1)
        Pnak = Pnak * 0.9;
    else if(celltype == 2)
        Pnak = Pnak * 0.7;

    *INaK = Pnak * (zna * JnakNa + zk * JnakK);
}

//Jrel
__inline__ __device__ void getJrel_ORd2011(real *Jrel, real *dJrelnp, real *dJrelp, real Jrelnp, real Jrelp, real ICaL, real cass, real cajsr, real fJrelp, int celltype, real Jrel_Multiplier) {

    real jsrMidpoint = 1.7f;

    real bt = 4.75f;
    real a_rel = 0.5f * bt;
    real Jrel_inf = a_rel * (-ICaL) / (1.0f + powf((jsrMidpoint / cajsr), 8.0f));

    if (celltype == 2)
        Jrel_inf = Jrel_inf * 1.7f;

    real tau_rel = bt / (1.0f + 0.0123f / cajsr);

    if (tau_rel < 0.001) tau_rel = 0.001f;

    *dJrelnp = (Jrel_inf - Jrelnp) / tau_rel;
    real btp = 1.25f * bt;
    real a_relp = 0.5f * btp;
    real Jrel_infp = a_relp * (-ICaL) / (1.0f + powf((jsrMidpoint / cajsr), 8.0f));

    if (celltype == 2)
        Jrel_infp = Jrel_infp * 1.7f;

    real tau_relp = btp / (1.0f + 0.0123f / cajsr);

    if (tau_relp < 0.001f)
        tau_relp = 0.001f;

    *dJrelp = (Jrel_infp - Jrelp) / tau_relp;

    *Jrel = Jrel_Multiplier * 1.5378f * ((1.0f - fJrelp) * Jrelnp + fJrelp * Jrelp);
}

//Jup
__inline__ __device__ void getJup_ORd2011(real *Jup, real *Jleak, real cai, real cansr, real fJupp, int celltype, real Jup_Multiplier) {
    // calculate serca pump, ca uptake flux
    real Jupnp = Jup_Multiplier * 0.005425 * cai / (cai + 0.00092);
    real Jupp = Jup_Multiplier * 2.75 * 0.005425 * cai / (cai + 0.00092 - 0.00017);
    if (celltype == 1) {
        Jupnp = Jupnp * 1.3f;
        Jupp = Jupp * 1.3f;
    }

    *Jleak = Jup_Multiplier * 0.0048825 * cansr / 15.0;
    *Jup = (1.0 - fJupp) * Jupnp + fJupp * Jupp - *Jleak;
}


