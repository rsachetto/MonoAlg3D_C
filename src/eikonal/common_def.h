//
// CUDA implementation of FIM (Fast Iterative Method) for Eikonal equations
//
// Copyright (c) Won-Ki Jeong (wkjeong@unist.ac.kr)
//
// 2016. 2. 4
//

//
// Common to entire project
//

#ifndef __COMMON_DEF_H__
#define __COMMON_DEF_H__

#include <stdbool.h>
#include "../common_types/common_types.h"

//
// common definition for Eikonal solvers
//
#ifndef INF
#define INF 1e20//FLT_MAX //
#endif

#define BLOCK_LENGTH 4

//
// itk image volume definition for 3D anisotropic eikonal solvers
//
typedef unsigned int uint;
typedef unsigned char uchar;

struct CUDA_MEM_STRUCTURE {
  // volsize/blksize : # of pixel in volume/block
  // blknum : # of block
  // blklength : # of pixel in one dimemsion of block
  uint nActiveBlock, blknum, volsize, blksize;
  int xdim, ydim, zdim, nIter, blklength; // new new x,y,z dim to aligh power of 4

  // host memory
  uint *h_list;
  bool *h_listVol, *h_listed;

  // device memory
  uint *d_list;
  real *d_spd;
  bool *d_mask, *d_listVol, *d_con;
  
  real *h_sol;//h_speedtable[256];
  real *d_sol, *t_sol; 

  // GroupOrder
  int* blockOrder;
  int K;
};

typedef struct CUDA_MEM_STRUCTURE CUDAMEMSTRUCT;

#endif
