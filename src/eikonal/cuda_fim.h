//
// CUDA implementation of FIM (Fast Iterative Method) for Eikonal equations
//
// Copyright (c) Won-Ki Jeong (wkjeong@unist.ac.kr)
//
// 2016. 2. 4
//
#ifndef __CUDA_FIM_H__
#define __CUDA_FIM_H__

#include "common_def.h"

#ifdef __cplusplus
extern "C" {
#endif
void run_eikonal_solver_simple(CUDAMEMSTRUCT *cmem, bool verbose);
#ifdef __cplusplus
}
#endif

#endif
