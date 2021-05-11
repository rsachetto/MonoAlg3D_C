#include <stddef.h>
#include <stdint.h>

#include "stewart_aslanidi_noble_2009.h"

__global__ void kernel_set_model_inital_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        *((real * )((char *) sv + pitch * 0) + threadID) = -69.1370441635924;
        *((real * )((char *) sv + pitch * 1) + threadID) = 136.781894160227;
        *((real * )((char *) sv + pitch * 2) + threadID) = 8.80420286531673;
        *((real * )((char *) sv + pitch * 3) + threadID) = 0.000101878186157052;
        *((real * )((char *) sv + pitch * 4) + threadID) = 0.0457562667986602;
        *((real * )((char *) sv + pitch * 5) + threadID) = 0.00550281999719088;
        *((real * )((char *) sv + pitch * 6) + threadID) = 0.313213286437995;
        *((real * )((char *) sv + pitch * 7) + threadID) = 0.00953708522974789;
        *((real * )((char *) sv + pitch * 8) + threadID) = 0.0417391656294997;
        *((real * )((char *) sv + pitch * 9) + threadID) = 0.190678733735145;
        *((real * )((char *) sv + pitch * 10) + threadID) = 0.238219836154029;
        *((real * )((char *) sv + pitch * 11) + threadID) = 0.000446818714055411;
        *((real * )((char *) sv + pitch * 12) + threadID) = 0.000287906256206415;
        *((real * )((char *) sv + pitch * 13) + threadID) = 0.989328560287987;
        *((real * )((char *) sv + pitch * 14) + threadID) = 0.995474890442185;
        *((real * )((char *) sv + pitch * 15) + threadID) = 0.999955429598213;
        *((real * )((char *) sv + pitch * 16) + threadID) = 0.96386101799501;
        *((real * )((char *) sv + pitch * 17) + threadID) = 0.00103618091196912;
        *((real * )((char *) sv + pitch * 18) + threadID) = 3.10836886659417;
        *((real * )((char *) sv + pitch * 19) + threadID) = 0.991580051907845;

        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * NEQ) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * (NEQ + 1)) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * (NEQ + 2)) + threadID) = 0.0;    // previous dt
        }
    }
}

inline __device__ void RHS_gpu(real *sv_, real *rDY_, real stim_current, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real STATES[NEQ];
    if (use_adpt_dt)
    {
        for (uint32_t i = 0; i < NEQ; i++)
            STATES[i] = sv_[i];
    }
    else
    {
        for (uint32_t i = 0; i < NEQ; i++)
            STATES[i] = *((real *)((char *)sv_ + pitch * i) + threadID_);
    }

    #include "stewart_aslanidi_noble_2009_common.inc"
}


//Include the default solver used by all current models.
#include "../default_solvers.cu"
