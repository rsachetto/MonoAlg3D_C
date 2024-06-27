//
// Created by sachetto on 06/06/24.
// Based on https://github.com/SCIInstitute/StructuredEikonal

#include <assert.h>
#include <cuda_runtime.h>
#include "eikonal_solver.h"
#include "../gpu_utils/gpu_utils.h"
#include "../logger/logger.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../3dparty/stb_ds.h"
#include "cuda_fim.h"

struct eikonal_solver * new_eikonal_solver(bool verbose) {

    struct eikonal_solver *solver = MALLOC_ONE_TYPE(struct eikonal_solver);
    solver->verbose = verbose;
    solver->cuda_mem_created = false;
    solver->width = 256;
    solver->height = 256;
    solver->depth = 256;
    solver->iters_per_block = 10;
    solver->solver_type = 0;
    solver->speeds = NULL;
    solver->seeds = NULL;
    solver->mask = NULL;

    return solver;
}

void free_eikonal_solver(struct eikonal_solver *solver) {

    FREE_3D_ARRAY(solver->speeds, solver->width, solver->height);
    FREE_3D_ARRAY(solver->answer, solver->width, solver->height);
    FREE_3D_ARRAY(solver->mask, solver->width, solver->height);

    if(solver->cuda_mem_created) {
        free(solver->memory_struct.h_sol);
        free(solver->memory_struct.h_list);
        free(solver->memory_struct.h_listed);
        free(solver->memory_struct.h_listVol);
        free(solver->memory_struct.blockOrder);
        check_cuda_error( cudaFree(solver->memory_struct.d_spd) );
        check_cuda_error( cudaFree(solver->memory_struct.d_sol) );
        check_cuda_error( cudaFree(solver->memory_struct.t_sol) );  // temp solution for ping-pong
        check_cuda_error( cudaFree(solver->memory_struct.d_con) );  // convergence volume
        check_cuda_error( cudaFree(solver->memory_struct.d_list) );
        check_cuda_error( cudaFree(solver->memory_struct.d_listVol) );
        check_cuda_error( cudaFree(solver->memory_struct.d_mask) );
    }

    free(solver);
}

void check_cuda_memory(struct eikonal_solver *solver) {
    if (solver->verbose) {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);

        printf("Total Memory : %lld MB\n", total_mem / (1024LL * 1024LL));
        printf("Free Memory  : %lld MB\n", free_mem / (1024LL * 1024LL));
        printf("--\n");
    }
}

void init_cuda_mem(struct eikonal_solver *solver) {

    if(solver->width <= 0 || solver->height <= 0 || solver->depth <= 0) {
        log_error_and_exit("Volume dimension cannot be zero");
    }

    check_cuda_memory(solver);

    // 1. Create /initialize GPU memory
    size_t nx, ny, nz;

    nx = solver->width + (BLOCK_LENGTH-solver->width%BLOCK_LENGTH)%BLOCK_LENGTH;
    ny = solver->height + (BLOCK_LENGTH-solver->height%BLOCK_LENGTH)%BLOCK_LENGTH;
    nz = solver->depth + (BLOCK_LENGTH-solver->depth%BLOCK_LENGTH)%BLOCK_LENGTH;

    if (solver->verbose) {
        printf("%ld %ld %ld \n",nx,ny,nz);
    }

    size_t volSize = nx*ny*nz;
    size_t blkSize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH;

    size_t nBlkX = nx / BLOCK_LENGTH;
    size_t nBlkY = ny / BLOCK_LENGTH;
    size_t nBlkZ = nz / BLOCK_LENGTH;
    size_t blockNum = nBlkX*nBlkY*nBlkZ;

    solver->memory_struct.xdim = (int) nx;
    solver->memory_struct.ydim = (int) ny;
    solver->memory_struct.zdim = (int) nz;
    solver->memory_struct.volsize = (uint) volSize;
    solver->memory_struct.blksize = (uint) blkSize;
    solver->memory_struct.blklength = BLOCK_LENGTH;
    solver->memory_struct.blknum = (uint) blockNum;
    solver->memory_struct.nIter = (int) solver->iters_per_block; // iter per block

    if(solver->cuda_mem_created) // delete previous memory
    {
        free(solver->memory_struct.h_sol);
        free(solver->memory_struct.h_list);
        free(solver->memory_struct.h_listed);
        free(solver->memory_struct.h_listVol);
        free(solver->memory_struct.blockOrder);
        check_cuda_error( cudaFree(solver->memory_struct.d_spd) );
        check_cuda_error( cudaFree(solver->memory_struct.d_sol) );
        check_cuda_error( cudaFree(solver->memory_struct.t_sol) );  // temp solution for ping-pong
        check_cuda_error( cudaFree(solver->memory_struct.d_con) );  // convergence volume
        check_cuda_error( cudaFree(solver->memory_struct.d_list) );
        check_cuda_error( cudaFree(solver->memory_struct.d_listVol) );
        check_cuda_error( cudaFree(solver->memory_struct.d_mask) );
    }
    solver->cuda_mem_created = true;

    solver->memory_struct.h_sol = (real*) malloc(volSize*sizeof(real)); // initial solution
    solver->memory_struct.h_list = (uint*) malloc(blockNum*sizeof(uint)); // linear list contains active block indices
    solver->memory_struct.h_listed = (bool*) malloc(blockNum*sizeof(bool));  // whether block is added to the list
    solver->memory_struct.h_listVol = (bool*) malloc(blockNum*sizeof(bool)); // volume list shows active/nonactive of corresponding block
    solver->memory_struct.blockOrder = (int*) malloc(blockNum*sizeof(int));

    check_cuda_memory(solver);

    //
    // create host/device memory using CUDA mem functions
    //
    check_cuda_error( cudaMalloc((void**)&(solver->memory_struct.d_spd), volSize*sizeof(real)) );
    check_cuda_memory(solver);

    check_cuda_error( cudaMalloc((void**)&(solver->memory_struct.d_sol), volSize*sizeof(real)) );
    check_cuda_memory(solver);

    check_cuda_error( cudaMalloc((void**)&(solver->memory_struct.t_sol), volSize*sizeof(real)) );  // temp solution for ping-pong
    check_cuda_memory(solver);

    check_cuda_error( cudaMalloc((void**)&(solver->memory_struct.d_con), volSize*sizeof(bool))  );  // convergence volume
    check_cuda_memory(solver);

    check_cuda_error( cudaMalloc((void**)&(solver->memory_struct.d_list), blockNum*sizeof(uint)) );
    check_cuda_memory(solver);

    check_cuda_error( cudaMalloc((void**)&(solver->memory_struct.d_listVol), blockNum*sizeof(bool)) );
    check_cuda_memory(solver);

    check_cuda_error( cudaMalloc((void**)&(solver->memory_struct.d_mask), volSize*sizeof(bool)) );
    check_cuda_memory(solver);
}


void set_attribute_mask(struct eikonal_solver *solver) {

    assert(solver);
    assert(solver->speeds);
    assert(solver->mask);

    uint volSize = solver->memory_struct.volsize;

    int nx, ny, nz, blklength;

    nx = solver->memory_struct.xdim;
    ny = solver->memory_struct.ydim;
    nz = solver->memory_struct.zdim;
    blklength = solver->memory_struct.blklength;

    // create host memory
    real *h_spd  = MALLOC_ARRAY_OF_TYPE(real, volSize); // byte speed, host
    bool  *h_mask = MALLOC_ARRAY_OF_TYPE(bool, volSize);

    // copy input volume to host memory
    // make each block to be stored contiguously in 1D memory space
    uint idx = 0;
    for(int zStr = 0; zStr < nz; zStr += blklength) {
        for(int yStr = 0; yStr < ny; yStr += blklength) {
            for(int xStr = 0; xStr < nx; xStr += blklength) {
                // for each block
                for(int z = zStr; z < zStr + blklength; z++) {
                    for(int y = yStr; y < yStr + blklength; y++) {
                        for(int x = xStr; x < xStr+blklength; x++) {
                            h_spd[idx] = solver->speeds[x][y][z];
                            h_mask[idx] = solver->mask[x][y][z];
                            idx++;
                        }
                    }
                }
            }
        }
    }

    // initialize GPU memory with host memory
    check_cuda_error( cudaMemcpy(solver->memory_struct.d_spd, h_spd, volSize*sizeof(real), cudaMemcpyHostToDevice) );
    check_cuda_error( cudaMemcpy(solver->memory_struct.d_mask, h_mask, volSize*sizeof(bool), cudaMemcpyHostToDevice) );

    free(h_spd);
    free(h_mask);
}

static void init_attribute_mask(struct eikonal_solver *solver) {

    size_t size_x = solver->width;
    size_t size_y = solver->height;
    size_t size_z = solver->depth;

    ALLOCATE_3D_ARRAY(solver->mask, bool, size_x, size_y, size_z);

    size_t n_active = solver->num_active_cells;
    struct cell_node **active_cells = solver->active_cells;

    for(size_t i = 0; i < n_active; i++) {
        int x = (int) active_cells[i]->center.x;
        int y = (int) active_cells[i]->center.y;
        int z = (int) active_cells[i]->center.z;
        solver->mask[x][y][z] = true;
    }

}

static void initialization(struct eikonal_solver *solver) {
    assert(solver);
    check_cuda_memory(solver);
    init_cuda_mem(solver);
    init_attribute_mask(solver);
    set_attribute_mask(solver);
    check_cuda_memory(solver);
}

static void get_solution(struct eikonal_solver *solver) {
    // copy solution from GPU
    check_cuda_error( cudaMemcpy(solver->memory_struct.h_sol,
                                 solver->memory_struct.d_sol,
                                 solver->memory_struct.volsize*sizeof(real),
                                 cudaMemcpyDeviceToHost) );

    size_t size_x = solver->width;
    size_t size_y = solver->height;
    size_t size_z = solver->depth;

    ALLOCATE_3D_ARRAY(solver->answer, real, size_x, size_y, size_z);

    for(size_t blockID = 0; blockID < solver->memory_struct.blknum; blockID++) {
        size_t baseAddr = blockID * solver->memory_struct.blksize;
        size_t xgridlength = solver->memory_struct.xdim/BLOCK_LENGTH;
        size_t ygridlength = solver->memory_struct.ydim/BLOCK_LENGTH;
        // compute block index
        size_t bx = blockID%xgridlength;
        size_t tmpIdx = (blockID - bx)/xgridlength;
        size_t by = tmpIdx%ygridlength;
        size_t bz = (tmpIdx-by)/ygridlength;
        //translate back to real space
        for(int k = 0; k < BLOCK_LENGTH; k++) {
            for(int j = 0; j < BLOCK_LENGTH; j++) {
                for(int i = 0; i < BLOCK_LENGTH; i++) {
                    double d = solver->memory_struct.h_sol[baseAddr + k * BLOCK_LENGTH * BLOCK_LENGTH + j * BLOCK_LENGTH + i];
                    if ((i + bx * BLOCK_LENGTH) < solver->width &&
                        (j + by * BLOCK_LENGTH) < solver->height &&
                        (k + bz * BLOCK_LENGTH) < solver->depth) {
                        solver->answer[(i + bx * BLOCK_LENGTH)][(j + by * BLOCK_LENGTH)][k + bz * BLOCK_LENGTH] = d;
                    }
                }
            }
        }
    }

    for(int i = 0 ; i < solver->num_active_cells; i++) {
        solver->active_cells[i]->v = solver->answer[(int)solver->active_cells[i]->center.x][(int)solver->active_cells[i]->center.y][(int)solver->active_cells[i]->center.z];

        //Translating back to original space
        solver->active_cells[i]->center.x = solver->active_cells[i]->center.x * solver->active_cells[i]->discretization.x + solver->active_cells[i]->discretization.x/2.0;
        solver->active_cells[i]->center.y  = solver->active_cells[i]->center.y * solver->active_cells[i]->discretization.y + solver->active_cells[i]->discretization.y/2.0;
        solver->active_cells[i]->center.z = solver->active_cells[i]->center.z * solver->active_cells[i]->discretization.z + solver->active_cells[i]->discretization.z/2.0;

    }
}

void map_generator(struct eikonal_solver *solver) {

    double pi = 3.141592653589793238462643383;

    size_t size_x = solver->width;
    size_t size_y = solver->height;
    size_t size_z = solver->depth;

    ALLOCATE_3D_ARRAY(solver->speeds, real, size_x, size_y, size_z);

    switch(solver->solver_type) {
    case 0 :
        //Constant Speed Map
        for (int k = 0 ; k < solver->depth ; ++k) {
            for (int j = 0 ; j < solver->height; ++j) {
                for ( int i = 0 ; i < solver->width ; ++i) {
                    solver->speeds[i][j][k] = 1.0;
                }
            }
        }
        break;
    case 1 :
        //Sinusoid Speed Map
        for (int k = 0 ; k < solver->depth ; ++k) {
            for (int j = 0 ; j < solver->height; ++j) {
                for ( int i = 0 ; i < solver->width ; ++i) {
                    solver->speeds[i][j][k] =
                        (6.0 + 5.0 *(sin((i*pi)/(double)solver->width * 2.0))*
                             sin((j*pi)/(double)solver->height*2.0)*
                             sin((k*pi)/(double)solver->depth*2.0));
                }
            }
        }

        break;
    }
}

void use_seeds(struct eikonal_solver *solver) {

    if (solver->verbose) {
        printf("Loading seed volume...\n");
    }
    uint volSize, blockNum;
    int nx, ny, nz, blklength;

    nx = solver->memory_struct.xdim;
    ny = solver->memory_struct.ydim;
    nz = solver->memory_struct.zdim;
    volSize = solver->memory_struct.volsize;
    blklength = solver->memory_struct.blklength;
    blockNum = solver->memory_struct.blknum;

    // copy input volume to host memory
    // make each block to be stored contiguously in 1D memory space
    uint idx = 0;
    uint blk_idx = 0;
    uint list_idx = 0;
    uint nActiveBlock = 0;

    for(int zStr = 0; zStr < nz; zStr += blklength) {
        for(int yStr = 0; yStr < ny; yStr += blklength) {
            for(int xStr = 0; xStr < nx; xStr += blklength) {
                // for each block
                bool isSeedBlock = false;

                for(int z=zStr; z<zStr+blklength; z++) {
                    for(int y=yStr; y<yStr+blklength; y++) {
                        for(int x=xStr; x<xStr+blklength; x++) {
                            solver->memory_struct.h_sol[idx] = INF;
                            if (arrlen(solver->seeds) == 0) {
                                if (x == nx/2 && y == ny/2 && z == nz/2) {
                                    solver->memory_struct.h_sol[idx] = 0;
                                    isSeedBlock = true;
                                    if (solver->verbose) {
                                        printf("%d is Selected bt source \n",idx);
                                    }
                                }
                            } else {
                                for(size_t i = 0; i < arrlen(solver->seeds); i++) {
                                    if (solver->seeds[i][0] == x &&
                                        solver->seeds[i][1] == y &&
                                        solver->seeds[i][2] == z) {
                                        solver->memory_struct.h_sol[idx] = 0;
                                        isSeedBlock = true;
                                        if (solver->verbose) {
                                            printf("%d is Selected bt source \n",idx);
                                        }
                                    }
                                }
                            }
                            idx++;
                        }
                    }
                }
                ///////////////////////////////////////////////
                if(isSeedBlock) {
                    if (solver->verbose) {
                        printf("%d,%d,%d is Seed Block \n",zStr,yStr,xStr);
                    }
                    solver->memory_struct.h_listVol[blk_idx] = true;
                    solver->memory_struct.h_listed[blk_idx] = true;
                    solver->memory_struct.h_list[list_idx] = blk_idx;
                    list_idx++;
                    nActiveBlock++;
                } else {
                    solver->memory_struct.h_listVol[blk_idx] = false;
                    solver->memory_struct.h_listed[blk_idx] = false;
                }
                blk_idx++;
            }
        }
    }
    solver->memory_struct.nActiveBlock = nActiveBlock;
    // initialize GPU memory with host memory
    check_cuda_error( cudaMemcpy(solver->memory_struct.d_sol, solver->memory_struct.h_sol, volSize*sizeof(real), cudaMemcpyHostToDevice) );
    check_cuda_error( cudaMemcpy(solver->memory_struct.t_sol, solver->memory_struct.h_sol, volSize*sizeof(real), cudaMemcpyHostToDevice) );
    check_cuda_error( cudaMemcpy(solver->memory_struct.d_list, solver->memory_struct.h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice) );
    check_cuda_error( cudaMemcpy(solver->memory_struct.d_listVol, solver->memory_struct.h_listVol, blockNum*sizeof(bool), cudaMemcpyHostToDevice) );
    // initialize GPU memory with constant value
    check_cuda_error( cudaMemset(solver->memory_struct.d_con, 1, volSize*sizeof(bool)) );
}

void solve_eikonal(struct eikonal_solver *solver) {

    if (solver->speeds == NULL) {
        map_generator(solver);
    }

    solver->cuda_mem_created = false;
    initialization(solver);
    use_seeds(solver);
    run_eikonal_solver_simple(&solver->memory_struct, solver->verbose);
    get_solution(solver);

}
