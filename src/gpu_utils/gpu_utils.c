#include "gpu_utils.h"

#include <stdlib.h>
#include <stdio.h>

void cuda_assert(cudaError_t code, const char *file, int line)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPU Error!: %s %s %d\n", cudaGetErrorString(code), file, line);
        exit(code);
    }
}

