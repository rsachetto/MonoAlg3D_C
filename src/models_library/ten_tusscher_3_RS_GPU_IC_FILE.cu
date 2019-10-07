#include <stddef.h>
#include "../monodomain/constants.h"
#include "model_gpu_utils.h"
extern "C" {
    #include "../string/sds.h"
}

#include "../single_file_libraries/stb_ds.h"

#include "ten_tusscher_3_RS.h"


extern "C" bool get_ic_from_file(real *params, real *ICs) {

    string_array lines = read_lines("ICs_1000configs.txt");

    size_t num_lines = arrlen(lines);
    int tmp;

    bool found = false;

    for(size_t i = 0; i < num_lines; i++) {

        sds *splitted_line = sdssplit(lines[i], " | ", &tmp);

        int num_par;
        sds *parameters = NULL;

        if(splitted_line[2][0] == '1') {

            parameters = sdssplit(splitted_line[1], " ", &num_par);

            for(int p = 0; p < num_par; p++) {
                
                real value = strtof(parameters[p], NULL);

                if(params[p] != value)  { 
                    printf("%lf %lf\n", params[p], value);
                    break;
                }

                found = true;
            
            }

            if (found) { 
                
                int num_ICs;
                sds *ICs_string = sdssplit(splitted_line[4], " ", &num_ICs);
                
                real V = strtof(splitted_line[3], NULL);

               // printf("%d\n", num_ICs);
               // printf("%lf\n", V);
                ICs[0] = V;

                for(int j = 0; j < NEQ-1; j++) {
                    ICs[j+1] =  strtof(ICs_string[j], NULL);
                    //printf("%s - %g\n", ICs_string[j], ICs[j+1]);
                }

                sdsfreesplitres(splitted_line, tmp);  
                sdsfreesplitres(parameters, num_par);  
                break; 
            }

        }

        sdsfreesplitres(splitted_line, tmp);  
        if(parameters) sdsfreesplitres(parameters, num_par);  

    }

    for(size_t i = 0; i < num_lines; i++) {
        free(lines[i]);   
    }

    arrfree(lines);
    return found;
    
}

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    char *cell_type;
    #ifdef ENDO
        cell_type = strdup("ENDO");
    #endif

    #ifdef EPI
        cell_type = strdup("EPI");
    #endif

    #ifdef MCELL
        cell_type = strdup("MCELL");
    #endif

    print_to_stdout_and_file("Using ten Tusscher 3 %s GPU model\n", cell_type);

    free(cell_type);

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);

    check_cuda_error(cudaMallocPitch((void **) &(*sv), &pitch_h, size, (size_t )NEQ));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));    
    real *ICs_device = NULL;
    
    if(extra_data) {        
        size_t mem = NEQ*sizeof(real);         
        real *ICs = (real*) malloc(mem);
        
        if(get_ic_from_file((real*)extra_data, ICs)) {            
            check_cuda_error(cudaMalloc((void **)&ICs_device, mem));
            check_cuda_error(cudaMemcpy(ICs_device, ICs, mem, cudaMemcpyHostToDevice));        
            free(ICs);

        }
        else {
            print_to_stderr_and_file_and_exit("Combination not found: %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",  
                                               ((real*)extra_data)[0], ((real*)extra_data)[1],((real*)extra_data)[2],
                                               ((real*)extra_data)[3], ((real*)extra_data)[4], ((real*)extra_data)[5], 
                                               ((real*)extra_data)[6]);
            free(ICs);

            exit(0);
        }
    }

    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(*sv, ICs_device, num_volumes);

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

    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    // Default values for a healthy cell ///////////
    real atpi = 6.8f;
    real Ko = 5.4f;
    real Ki = 138.3f;
    real Vm_change = 0.0;
    real GNa_multiplicator = 1.0f;
    real GCaL_multiplicator = 1.0f;
    real INaCa_multiplicator = 1.0f;
    ////////////////////////////////////

    real *fibrosis_device;
    real *fibs = NULL;
    int num_extra_parameters = 7;
    size_t extra_parameters_size = num_extra_parameters*sizeof(real);

    real *extra_parameters_device;
    real fibs_size = num_cells_to_solve*sizeof(real);

    bool dealocate = false;

    if(extra_data) {
        fibs = ((real*)extra_data) + num_extra_parameters; //pointer
    }
    else {
        extra_data = malloc(extra_parameters_size);
        ((real*)extra_data)[0] = atpi;
        ((real*)extra_data)[1] = Ko;
        ((real*)extra_data)[2] = Ki;
        ((real*)extra_data)[3] = Vm_change;
        ((real*)extra_data)[4] = GNa_multiplicator;
        ((real*)extra_data)[5] = GCaL_multiplicator;
        ((real*)extra_data)[6] = INaCa_multiplicator;

        fibs = (real*)calloc(num_cells_to_solve, sizeof(real));

        dealocate = true;
    }

    check_cuda_error(cudaMalloc((void **) &extra_parameters_device, extra_parameters_size));
    check_cuda_error(cudaMemcpy(extra_parameters_device, extra_data, extra_parameters_size, cudaMemcpyHostToDevice));

    check_cuda_error(cudaMalloc((void **) &fibrosis_device, fibs_size));
    check_cuda_error(cudaMemcpy(fibrosis_device, fibs, fibs_size, cudaMemcpyHostToDevice));

    solve_gpu<<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, fibrosis_device, extra_parameters_device);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    check_cuda_error(cudaFree(fibrosis_device));
    check_cuda_error(cudaFree(extra_parameters_device));

    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));

    if(dealocate) {
        free(fibs);
        free(extra_data);
    }
}

__global__ void kernel_set_model_inital_conditions(real *sv, real *ICs, int num_volumes)
{
    // Thread ID
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if(threadID < num_volumes) {

        if(ICs == NULL) {
            *((real *) ((char *) sv + pitch * 0) + threadID) = INITIAL_V;   // V;       millivolt
            *((real *) ((char *) sv + pitch * 1) + threadID) = 0.0f; //M
            *((real *) ((char *) sv + pitch * 2) + threadID) = 0.75; //H
            *((real *) ((char *) sv + pitch * 3) + threadID) = 0.75; //J
            *((real *) ((char *) sv + pitch * 4) + threadID) = 0.0f; //Xr1
            *((real *) ((char *) sv + pitch * 5) + threadID) = 0.0f; //Xs
            *((real *) ((char *) sv + pitch * 6) + threadID) = 1.0; //S
            *((real *) ((char *) sv + pitch * 7) + threadID) = 1.0; //F
            *((real *) ((char *) sv + pitch * 8) + threadID) = 1.0; //F2
            *((real *) ((char *) sv + pitch * 9) + threadID) = 0.0; //D_INF
            *((real *) ((char *) sv + pitch * 10) + threadID) = 0.0; //R_INF
            *((real *) ((char *) sv + pitch * 11) + threadID) = 0.0; //Xr2_INF
        }
        else {
            *((real *) ((char *) sv + pitch * 0) + threadID) = ICs[0];   // V;       millivolt
            *((real *) ((char *) sv + pitch * 1) + threadID) = ICs[1]; //M
            *((real *) ((char *) sv + pitch * 2) + threadID) = ICs[2]; //H
            *((real *) ((char *) sv + pitch * 3) + threadID) = ICs[3]; //J
            *((real *) ((char *) sv + pitch * 4) + threadID) = ICs[4]; //Xr1
            *((real *) ((char *) sv + pitch * 5) + threadID) = ICs[5]; //Xs
            *((real *) ((char *) sv + pitch * 6) + threadID) = ICs[6]; //S
            *((real *) ((char *) sv + pitch * 7) + threadID) = ICs[7]; //F
            *((real *) ((char *) sv + pitch * 8) + threadID) = ICs[8]; //F2
            *((real *) ((char *) sv + pitch * 9) + threadID) = ICs[9]; //D_INF
            *((real *) ((char *) sv + pitch * 10) + threadID) = ICs[10]; //R_INF
            *((real *) ((char *) sv + pitch * 11) + threadID) = ICs[11]; //Xr2_INF
        }
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps, real *fibrosis, real *extra_parameters)
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

            RHS_gpu(sv, rDY, stim_currents[threadID], sv_id, dt, fibrosis[threadID], extra_parameters);

            *((real*)((char*)sv) + sv_id) = dt*rDY[0] + *((real*)((char*)sv) + sv_id);

            for(int i = 1; i < 12; i++) {
                *((real*)((char*)sv + pitch * i) + sv_id) = rDY[i];
            }

        }

    }
}


inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt, real fibrosis, real *extra_parameters) {

    //fibrosis = 0 means that the cell is fibrotic, 1 is not fibrotic. Anything between 0 and 1 means border zone
    //printf("%lf, %lf \n", extra_parameters[0], extra_parameters[1]);

    const real svolt = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real sm   = *((real*)((char*)sv_ + pitch * 1) + threadID_);
    const real sh   = *((real*)((char*)sv_ + pitch * 2) + threadID_);
    const real sj   = *((real*)((char*)sv_ + pitch * 3) + threadID_);
    const real sxr1 = *((real*)((char*)sv_ + pitch * 4) + threadID_);
    const real sxs  = *((real*)((char*)sv_ + pitch * 5) + threadID_);
    const real ss   = *((real*)((char*)sv_ + pitch * 6) + threadID_);
    const real sf  = *((real*)((char*)sv_ + pitch * 7) + threadID_);
    const real sf2  = *((real*)((char*)sv_ + pitch * 8) + threadID_);
    const real D_INF  = *((real*)((char*)sv_ + pitch * 9) + threadID_);
    const real R_INF  = *((real*)((char*)sv_ + pitch * 10) + threadID_);
    const real Xr2_INF  = *((real*)((char*)sv_ + pitch * 11) + threadID_);

    //FUCK YOU NVCC
    #include "ten_tusscher_3_RS_common.inc"


}