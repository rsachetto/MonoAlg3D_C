#include "../../gpu_utils/gpu_utils.h"
#include "noble_1962.h"
#include <stddef.h>
#include <stdint.h>

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    log_info("Using noble_1962 GPU model\n");

    uint32_t num_volumes = solver->original_num_cells;

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    // allocates a 2d contigous array
    check_cuda_error(cudaMallocPitch((void **) &(solver->sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(solver->sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
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
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));


    //the array cells to solve is passed when we are using and adaptive mesh
    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }
    solve_gpu <<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        // Default values
        //*((real * )((char *) sv + pitch * 0) + threadID) = -75.5344986658f; //V millivolt 
        //*((real * )((char *) sv + pitch * 1) + threadID) = 0.060546727200f;   //m dimensionless 
        //*((real * )((char *) sv + pitch * 2) + threadID) = 0.725900135500f;   //h millivolt 
        //*((real * )((char *) sv + pitch * 3) + threadID) = 0.470923970800f;   //n dimensionless 

        // BCL = 300ms
        *((real * )((char *) sv + pitch * 0) + threadID) = -81.1893;    // V millivolt 
        *((real * )((char *) sv + pitch * 1) + threadID) = 0.0443563;    // m dimensionless
        *((real * )((char *) sv + pitch * 2) + threadID) = 0.851652;    // h dimensionless
        *((real * )((char *) sv + pitch * 3) + threadID) = 0.58291;    // n dimensionless
         
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
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

        real rDY[NEQ];

        for (int n = 0; n < num_steps; ++n) {

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id);

            for(int i = 0; i < NEQ; i++) {
                *((real *) ((char *) sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *) ((char *) sv + pitch * i) + sv_id);
            }            

        }

    }
}

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_) {

    //State variables
    const real V_old_ =  *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real m_old_ =  *((real*)((char*)sv_ + pitch * 1) + threadID_);
    const real h_old_ =  *((real*)((char*)sv_ + pitch * 2) + threadID_);
    const real n_old_ =  *((real*)((char*)sv_ + pitch * 3) + threadID_);

    //Parameters
    //const real Cm = 12.00000000000000000e+00f;             // (microF)
    //const real g_na_max = 400000.00000000000000000e+00f;   // (microS)
    //const real E_na = 40.00000000000000000e+00f;           // (millivolt)
    //const real g_L = 75.00000000000000000e+00f;            // (microS)
    //const real E_L = -60.00000000000000000e+00f;           // (millivolt)
    const real Cm = 12.0;             // (microF)
    const real g_na_max = 400.0;   // (microS)
    const real E_na = 40.0;           // (millivolt)
    const real g_L = 0.075;            // (microS)
    const real E_L = -60.0;           // (millivolt)


    real calc_I_stim = stim_current;

    // Algebraics
    //real g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    //real alpha_m = ( 100.000*(- V_old_ - 48.0000))/(exp((- V_old_ - 48.0000)/15.0000) - 1.00000);
    //real alpha_h =  170.000*exp((- V_old_ - 90.0000)/20.0000);
    //real alpha_n = ( 0.100000*(- V_old_ - 50.0000))/(exp((- V_old_ - 50.0000)/10.0000) - 1.00000);
    //real i_na =  (g_na+140.000)*(V_old_ - E_na);
    //real i_na_no_oscilation = (g_na+122.500)*(V_old_ - E_na);
    //real beta_m = ( 120.000*(V_old_+8.00000))/(exp((V_old_+8.00000)/5.00000) - 1.00000);
    //real beta_h = 1000.00/(1.00000+exp((- V_old_ - 42.0000)/10.0000));
    //real beta_n =  2.00000*exp((- V_old_ - 90.0000)/80.0000);
    //real g_K1 =  1200.00*exp((- V_old_ - 90.0000)/50.0000)+ 15.0000*exp((V_old_+90.0000)/60.0000);
    //real g_K2 =  1200.00*pow(n_old_, 4.00000);
    //real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    //real i_leak =  g_L*(V_old_ - E_L);

    real g_na =  pow(m_old_, 3.00000)*h_old_*g_na_max;
    real alpha_h = ((1.7e-01*exp((((-V_old_)-9.0e+01)/2.0e+01))));
    real alpha_m = (((1.0e-01*((-V_old_)-4.8e+01))/(exp((((-V_old_)-4.8e+01)/1.5e+01))-1.0e+00)));
    real alpha_n = (((1.0e-04*((-V_old_)-5.0e+01))/(exp((((-V_old_)-5.0e+01)/1.0e+01))-1.0e+00)));
    real i_na = (g_na+1.4e-01)*(V_old_ - E_na);
    //real i_na_no_oscilation = (g_na+1.2e-01)*(V_old_ - E_na);
    double beta_m = (((1.2e-01*(V_old_+8.0e+00))/(exp(((V_old_+8.0e+00)/5.0e+00))-1.0e+00)));
    double beta_h = ((1.0/(1.0e+00+exp((((-V_old_)-4.2e+01)/1.0e+01)))));
    double beta_n = ((2.0e-03*exp((((-V_old_)-9.0e+01)/8.0e+01))));
    real g_K1 = 1.2*exp((((-V_old_)-9.0e+01)/5.0e+01)) + (1.5e-02*exp(((V_old_+9.0e+01)/6.0e+01)));
    real g_K2 = 1.2*pow(n_old_,4.0e+00);
    real i_k =  (g_K1+g_K2)*(V_old_+100.000);
    real i_leak =  g_L*(V_old_ - E_L);

    // Rates
    //rDY_[0] = (- (i_na + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm) * 1.0E-03;
    //rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  beta_m*m_old_) * 1.0E-03;
    //rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  beta_h*h_old_) * 1.0E-03;
    //rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  beta_n*n_old_) * 1.0E-03;

    rDY_[0] = (- (i_na + i_k + i_leak + calc_I_stim)/Cm);
    //rDY_[0] = (- (i_na_no_oscilation + i_k + i_leak + calc_I_stim)/Cm);
    rDY_[1] =  (alpha_m*(1.00000 - m_old_) -  beta_m*m_old_);
    rDY_[2] =  (alpha_h*(1.00000 - h_old_) -  beta_h*h_old_);
    rDY_[3] =  (alpha_n*(1.00000 - n_old_) -  beta_n*n_old_);

}

