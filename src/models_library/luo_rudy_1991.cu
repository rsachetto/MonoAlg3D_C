#include "luo_rudy_1991.h"
#include <stddef.h>
#include <stdint.h>
#include "model_gpu_utils.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    print_to_stdout_and_file("Using luo_rudy_1991 GPU model\n");

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));


    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}

extern "C" SOLVE_MODEL_ODES_GPU(solve_model_odes_gpu) {

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
    solve_gpu <<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

}

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

         *((real * )((char *) sv + pitch * 0) + threadID) = -84.380111f; //V millivolt 
         *((real * )((char *) sv + pitch * 1) + threadID) = 0.001713f; //m dimensionless 
         *((real * )((char *) sv + pitch * 2) + threadID) = 0.982661f; //h dimensionless 
         *((real * )((char *) sv + pitch * 3) + threadID) = 0.989108f; //j dimensionless 
         *((real * )((char *) sv + pitch * 4) + threadID) = 0.003021f; //d dimensionless 
         *((real * )((char *) sv + pitch * 5) + threadID) = 0.999968f; //f dimensionless 
         *((real * )((char *) sv + pitch * 6) + threadID) = 0.041760f; //X dimensionless 
         *((real * )((char *) sv + pitch * 7) + threadID) = 0.000179f; //Cai millimolar 
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


#define IFNUMBER_0(name)if(((time_new>=stim_start)&&(time_new<=stim_end)&&(((time_new-stim_start)-(floor(((time_new-stim_start)/stim_period))*stim_period))<=stim_duration))) { (name) = stim_amplitude;    }  else{ (name) = 0.000000000000000e+00f;    }
#define IFNUMBER_1(name)if((V_old_<(-4.000000000000000e+01f))) { (name) = (1.350000000000000e-01f*expf(((8.000000000000000e+01f+V_old_)/(-6.800000000000000e+00f))));    }  else{ (name) = 0.000000000000000e+00f;    }
#define IFNUMBER_2(name)if((V_old_<(-4.000000000000000e+01f))) { (name) = ((3.560000000000000e+00f*expf((7.900000000000000e-02f*V_old_)))+(3.100000000000000e+05f*expf((3.500000000000000e-01f*V_old_))));    }  else{ (name) = (1.000000000000000e+00f/(1.300000000000000e-01f*(1.000000000000000e+00f+expf(((V_old_+1.066000000000000e+01f)/(-1.110000000000000e+01f))))));    }
#define IFNUMBER_3(name)if((V_old_<(-4.000000000000000e+01f))) { (name) = (((((-1.271400000000000e+05f)*expf((2.444000000000000e-01f*V_old_)))-(3.474000000000000e-05f*expf(((-4.391000000000000e-02f)*V_old_))))*(V_old_+3.778000000000000e+01f))/(1.000000000000000e+00f+expf((3.110000000000000e-01f*(V_old_+7.923000000000000e+01f)))));    }  else{ (name) = 0.000000000000000e+00f;    }
#define IFNUMBER_4(name)if((V_old_<(-4.000000000000000e+01f))) { (name) = ((1.212000000000000e-01f*expf(((-1.052000000000000e-02f)*V_old_)))/(1.000000000000000e+00f+expf(((-1.378000000000000e-01f)*(V_old_+4.014000000000000e+01f)))));    }  else{ (name) = ((3.000000000000000e-01f*expf(((-2.535000000000000e-07f)*V_old_)))/(1.000000000000000e+00f+expf(((-1.000000000000000e-01f)*(V_old_+3.200000000000000e+01f)))));    }
#define IFNUMBER_5(name)if((V_old_>(-1.000000000000000e+02f))) { (name) = ((2.837000000000000e+00f*(expf((4.000000000000000e-02f*(V_old_+7.700000000000000e+01f)))-1.000000000000000e+00f))/((V_old_+7.700000000000000e+01f)*expf((4.000000000000000e-02f*(V_old_+3.500000000000000e+01f)))));    }  else{ (name) = 1.000000000000000e+00f;    }

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_) {

    //State variables
    const real V_old_ =  *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real m_old_ =  *((real*)((char*)sv_ + pitch * 1) + threadID_);
    const real h_old_ =  *((real*)((char*)sv_ + pitch * 2) + threadID_);
    const real j_old_ =  *((real*)((char*)sv_ + pitch * 3) + threadID_);
    const real d_old_ =  *((real*)((char*)sv_ + pitch * 4) + threadID_);
    const real f_old_ =  *((real*)((char*)sv_ + pitch * 5) + threadID_);
    const real X_old_ =  *((real*)((char*)sv_ + pitch * 6) + threadID_);
    const real Cai_old_ =  *((real*)((char*)sv_ + pitch * 7) + threadID_);

    //Parameters
    const real C = 1.000000000000000e+00f;
    const real R = 8.314000000000000e+03f;
    const real T = 3.100000000000000e+02f;
    const real F = 9.648460000000001e+04f;
    const real Nao = 1.400000000000000e+02f;
    const real Nai = 1.800000000000000e+01f;
    const real g_Na = 2.300000000000000e+01f;
    const real Ko = 5.400000000000000e+00f;
    const real PR_NaK = 1.833000000000000e-02f;
    const real Ki = 1.450000000000000e+02f;
    const real g_Kp = 1.830000000000000e-02f;
    const real g_b = 3.921000000000000e-02f;
    const real E_b = -5.987000000000000e+01f;

    real calc_I_stim = stim_current;
    real calc_E_Na = (((R*T)/F)*logf((Nao/Nai)));	//2
    real calc_alpha_m = ((3.200000000000000e-01f*(V_old_+4.713000000000000e+01f))/(1.000000000000000e+00f-expf(((-1.000000000000000e-01f)*(V_old_+4.713000000000000e+01f)))));	//4
    real calc_beta_m = (8.000000000000000e-02f*expf(((-V_old_)/1.100000000000000e+01f)));	//5
    real calc_alpha_h = 0.0f;
    IFNUMBER_1(calc_alpha_h);	//7
    real calc_beta_h = 0.0f;
    IFNUMBER_2(calc_beta_h);	//8
    real calc_alpha_j = 0.0f;
    IFNUMBER_3(calc_alpha_j);	//10
    real calc_beta_j = 0.0f;
    IFNUMBER_4(calc_beta_j);	//11
    real calc_E_si = (7.700000000000000e+00f-(1.302870000000000e+01f*logf((Cai_old_/1.000000000000000e+00f))));	//13
    real calc_alpha_d = ((9.500000000000000e-02f*expf(((-1.000000000000000e-02f)*(V_old_-5.000000000000000e+00f))))/(1.000000000000000e+00f+expf(((-7.199999999999999e-02f)*(V_old_-5.000000000000000e+00f)))));	//15
    real calc_beta_d = ((7.000000000000001e-02f*expf(((-1.700000000000000e-02f)*(V_old_+4.400000000000000e+01f))))/(1.000000000000000e+00f+expf((5.000000000000000e-02f*(V_old_+4.400000000000000e+01f)))));	//16
    real calc_alpha_f = ((1.200000000000000e-02f*expf(((-8.000000000000000e-03f)*(V_old_+2.800000000000000e+01f))))/(1.000000000000000e+00f+expf((1.500000000000000e-01f*(V_old_+2.800000000000000e+01f)))));	//18
    real calc_beta_f = ((6.500000000000000e-03f*expf(((-2.000000000000000e-02f)*(V_old_+3.000000000000000e+01f))))/(1.000000000000000e+00f+expf(((-2.000000000000000e-01f)*(V_old_+3.000000000000000e+01f)))));	//19
    real calc_g_K = (2.820000000000000e-01f*powf((Ko/5.400000000000000e+00f),1.0/2.0));	//21
    real calc_E_K = (((R*T)/F)*logf(((Ko+(PR_NaK*Nao))/(Ki+(PR_NaK*Nai)))));	//22
    real calc_alpha_X = ((5.000000000000000e-04f*expf((8.300000000000000e-02f*(V_old_+5.000000000000000e+01f))))/(1.000000000000000e+00f+expf((5.700000000000000e-02f*(V_old_+5.000000000000000e+01f)))));	//24
    real calc_beta_X = ((1.300000000000000e-03f*expf(((-6.000000000000000e-02f)*(V_old_+2.000000000000000e+01f))))/(1.000000000000000e+00f+expf(((-4.000000000000000e-02f)*(V_old_+2.000000000000000e+01f)))));	//25
    real calc_Xi = 0.0f;
    IFNUMBER_5(calc_Xi);	//27
    real calc_g_K1 = (6.047000000000000e-01f*powf((Ko/5.400000000000000e+00f),1.0/2.0));	//28
    real calc_E_K1 = (((R*T)/F)*logf((Ko/Ki)));	//29
    real calc_Kp = (1.000000000000000e+00f/(1.000000000000000e+00f+expf(((7.488000000000000e+00f-V_old_)/5.980000000000000e+00f))));	//35
    real calc_i_b = (g_b*(V_old_-E_b));	//37
    real calc_i_Na = (g_Na*powf(m_old_,3.000000000000000e+00f)*h_old_*j_old_*(V_old_-calc_E_Na));	//3
    real calc_i_si = (9.000000000000000e-02f*d_old_*f_old_*(V_old_-calc_E_si));	//14
    real calc_alpha_K1 = (1.020000000000000e+00f/(1.000000000000000e+00f+expf((2.385000000000000e-01f*((V_old_-calc_E_K1)-5.921500000000000e+01f)))));	//31
    real calc_beta_K1 = (((4.912400000000000e-01f*expf((8.032000000000000e-02f*((V_old_+5.476000000000000e+00f)-calc_E_K1))))+(1.000000000000000e+00f*expf((6.175000000000000e-02f*(V_old_-(calc_E_K1+5.943099999999999e+02f))))))/(1.000000000000000e+00f+expf(((-5.143000000000000e-01f)*((V_old_-calc_E_K1)+4.753000000000000e+00f)))));	//32
    real calc_E_Kp = calc_E_K1;	//34
    real calc_i_K = (calc_g_K*X_old_*calc_Xi*(V_old_-calc_E_K));	//23
    real calc_K1_infinity = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));	//33
    real calc_i_Kp = (g_Kp*calc_Kp*(V_old_-calc_E_Kp));	//36
    real calc_i_K1 = (calc_g_K1*calc_K1_infinity*(V_old_-calc_E_K1));	//30

    rDY_[0] = (((-1.000000000000000e+00f)/C)*(calc_I_stim+calc_i_Na+calc_i_si+calc_i_K+calc_i_K1+calc_i_Kp+calc_i_b));
    rDY_[1] = ((calc_alpha_m*(1.000000000000000e+00f-m_old_))-(calc_beta_m*m_old_));
    rDY_[2] = ((calc_alpha_h*(1.000000000000000e+00f-h_old_))-(calc_beta_h*h_old_));
    rDY_[3] = ((calc_alpha_j*(1.000000000000000e+00f-j_old_))-(calc_beta_j*j_old_));
    rDY_[4] = ((calc_alpha_d*(1.000000000000000e+00f-d_old_))-(calc_beta_d*d_old_));
    rDY_[5] = ((calc_alpha_f*(1.000000000000000e+00f-f_old_))-(calc_beta_f*f_old_));
    rDY_[6] = ((calc_alpha_X*(1.000000000000000e+00f-X_old_))-(calc_beta_X*X_old_));
    rDY_[7] = ((((-1.000000000000000e-04f)/1.000000000000000e+00f)*calc_i_si)+(7.000000000000001e-02f*(1.000000000000000e-04f-Cai_old_)));

}

