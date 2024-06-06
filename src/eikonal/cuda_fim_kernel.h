//
// CUDA implementation of FIM (Fast Iterative Method) for Eikonal equations
//
// Copyright (c) Won-Ki Jeong (wkjeong@unist.ac.kr)
//
// 2016. 2. 4
//

#ifndef _cuda_fim_KERNEL_H_
#define _cuda_fim_KERNEL_H_

#include "common_def.h"

#define MEM(index) _mem[index]
#define SOL(i,j,k) _sol[i][j][k]
#define SPD(i,j,k) _spd[i][j][k]

__device__ real get_time_eikonal(real a, real b, real c, real s);
//
// F : Input speed (positive)
// if F =< 0, skip that pixel (masking out)
//
__global__ void run_solver(real* spd, bool* mask, const real *sol_in, 
  real *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim,
  int nIter, uint nActiveBlock);
//
// run_reduction
//
// con is pixelwise convergence. Do reduction on active tiles and write tile-wise
// convergence to listVol. The implementation assumes that the block size is 4x4x4.
//
__global__ void run_reduction(bool *con, bool *listVol, uint *list, uint nActiveBlock);
//
// if block is active block, copy values
// if block is neighbor, run solver once
//
__global__ void run_check_neighbor(real* spd, bool* mask, const real *sol_in, real *sol_out,
  bool *con, uint* list, int xdim, int ydim, int zdim,
  uint nActiveBlock, uint nTotalBlock);
#endif // #ifndef _cuda_fim_KERNEL_H_

