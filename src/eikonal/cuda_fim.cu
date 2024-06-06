//
// CUDA implementation of FIM (Fast Iterative Method) for Eikonal equations
//
// Copyright (c) Won-Ki Jeong (wkjeong@unist.ac.kr)
//
// 2016. 2. 4
//

#include <cstdio>
#include <assert.h>
#include "cuda_fim_kernel.h"
#include "cuda_fim.h"
#include "../gpu_utils/gpu_utils.h"
#include <math.h>

void run_eikonal_solver_simple(CUDAMEMSTRUCT *cmem, bool verbose) {
    
    int deviceID;
    cudaGetDevice(&deviceID);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, deviceID);
    if (verbose) {
        printf("Current device id : %d, name : %s\n", deviceID, deviceProp.name);
    }

    int xdim, ydim, zdim;
    xdim = cmem->xdim;
    ydim = cmem->ydim;
    zdim = cmem->zdim;

    // create volumes
    uint volSize = cmem->volsize;
    uint blockNum = cmem->blknum;

    if (verbose) {
        printf("# of total voxels : %d\n", volSize);
        printf("# of total blocks : %d\n", blockNum);
    }

    // h_ : host memory, d_ : device memory
    int nIter = cmem->nIter;
    uint nActiveBlock = cmem->nActiveBlock; // active list

    real*d_spd = cmem->d_spd;
    real *d_sol = cmem->d_sol;
    real *t_sol = cmem->t_sol;

    uint *d_list = cmem->d_list;
    bool *d_listVol = cmem->d_listVol;

    bool *d_con = cmem->d_con;
    bool *d_mask = cmem->d_mask;

    // copy so that original value should not be modified
    uint *h_list = (uint*) malloc(blockNum*sizeof(uint));
    bool *h_listed = (bool*) malloc(blockNum*sizeof(bool));
    bool *h_listVol = (bool*) malloc(blockNum*sizeof(bool));

    // initialization
    memcpy(h_list, cmem->h_list, blockNum*sizeof(uint));
    memcpy(h_listed, cmem->h_listed, blockNum*sizeof(bool));
    memcpy(h_listVol, cmem->h_listVol, blockNum*sizeof(bool));

    check_cuda_error( cudaMemcpy(cmem->d_list, cmem->h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice) );
    check_cuda_error( cudaMemcpy(cmem->d_listVol, cmem->h_listVol, blockNum*sizeof(bool), cudaMemcpyHostToDevice) );
    check_cuda_error( cudaMemcpy(cmem->d_sol, cmem->h_sol, volSize*sizeof(real), cudaMemcpyHostToDevice) );
    check_cuda_error( cudaMemcpy(cmem->t_sol, cmem->h_sol, volSize*sizeof(real), cudaMemcpyHostToDevice) );
    check_cuda_error( cudaMemset(cmem->d_con, 1, volSize*sizeof(bool)) );

    // set dimension of block and entire grid size
    dim3 dimBlock(BLOCK_LENGTH,BLOCK_LENGTH,BLOCK_LENGTH);
    dim3 dimEntireGrid(blockNum);
    dim3 dimGrid(nActiveBlock);

    int nTotalIter = 0;

    uint nTotalBlockProcessed = 0;

    // start solver
    while(nActiveBlock > 0)
    {
        assert(nActiveBlock < 4294967295);

        nTotalBlockProcessed += nActiveBlock;

        nTotalIter++;

        //
        // solve current blocks in the active lists
        //

        //            printf("# of active tiles : %u\n", nActiveBlock);
        if (verbose) {
            printf("# of active tiles : %u\n", nActiveBlock);
        }
        //////////////////////////////////////////////////////////////////
        // 1. run solver on current active tiles

        dimGrid.y = (unsigned int)floor(((double)nActiveBlock-1)/65535)+1;
        dimGrid.x = (unsigned int)ceil ((double)nActiveBlock/(double)dimGrid.y);

        if (verbose) {
            printf("Grid size : %d x %d\n", dimGrid.x, dimGrid.y);
        }

        check_cuda_error( cudaMemcpy(d_list, h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice) );

        run_solver<<< dimGrid, dimBlock >>>(d_spd, d_mask, d_sol, t_sol, d_con, d_list, xdim, ydim, zdim, nIter, nActiveBlock);

        check_cuda_error(cudaGetLastError());

        cudaDeviceSynchronize();


        //////////////////////////////////////////////////////////////////
        // 2. reduction (only active tiles)

        run_reduction<<< dimGrid, dim3(BLOCK_LENGTH,BLOCK_LENGTH,BLOCK_LENGTH/2) >>>(d_con, d_listVol, d_list, nActiveBlock);

        check_cuda_error(cudaGetLastError());

        cudaDeviceSynchronize();

        //////////////////////////////////////////////////////////////////
        // 3. check neighbor tiles of converged tile
        // Add any active block of neighbor of converged block is inserted
        // to the list

        check_cuda_error( cudaMemcpy(h_listVol, d_listVol, blockNum*sizeof(bool), cudaMemcpyDeviceToHost) );

        uint nOldActiveBlock = nActiveBlock;
        uint nBlkX = xdim/BLOCK_LENGTH;
        uint nBlkY = ydim/BLOCK_LENGTH;

        for(uint i=0; i<nOldActiveBlock; i++)
        {
            // check 6-neighbor of current active tile
            uint currBlkIdx = h_list[i];

            if(!h_listVol[currBlkIdx]) // not active : converged
            {
                uint nb[6];
                nb[0] = (currBlkIdx < nBlkX*nBlkY) ? currBlkIdx : (currBlkIdx - nBlkX*nBlkY);    //tp
                nb[1] = ((currBlkIdx + nBlkX*nBlkY) >= blockNum) ? currBlkIdx : (currBlkIdx + nBlkX*nBlkY); //bt
                nb[2] = (currBlkIdx < nBlkX) ? currBlkIdx : (currBlkIdx - nBlkX); //up
                nb[3] = ((currBlkIdx + nBlkX) >= blockNum) ? currBlkIdx : (currBlkIdx + nBlkX); //dn
                nb[4] = (currBlkIdx%nBlkX == 0) ? currBlkIdx : currBlkIdx-1; //lf
                nb[5] = ((currBlkIdx+1)%nBlkX == 0) ? currBlkIdx : currBlkIdx+1; //rt

                for(int nbIdx = 0; nbIdx < 6; nbIdx++)
                {
                    uint currIdx = nb[nbIdx];

                    if(!h_listed[currIdx])
                    {
                        h_listed[currIdx] = true;
                        h_list[nActiveBlock++] = currIdx;
                    }
                }
            }
        }
        cudaDeviceSynchronize();

        //////////////////////////////////////////////////////////////////
        // 4. run solver only once for neighbor blocks of converged block
        // current active list contains active blocks and neighbor blocks of
        // any converged blocks.
        //

        // update grid dimension because nActiveBlock is changed
        dimGrid.y = (unsigned int)floor(((double)nActiveBlock-1)/65535)+1;
        dimGrid.x = (unsigned int)ceil((double)nActiveBlock/(double)dimGrid.y);

        if (verbose) {
            printf("Grid size : %d x %d\n", dimGrid.x, dimGrid.y);
        }

        check_cuda_error(cudaMemcpy(d_list, h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice) );
        run_check_neighbor<<< dimGrid, dimBlock >>>(d_spd, d_mask, t_sol, d_sol, d_con, d_list, xdim, ydim, zdim, nOldActiveBlock, nActiveBlock);
        check_cuda_error(cudaGetLastError());
        cudaDeviceSynchronize();


        //////////////////////////////////////////////////////////////////
        // 5. reduction

        run_reduction<<< dimGrid, dim3(BLOCK_LENGTH,BLOCK_LENGTH,BLOCK_LENGTH/2) >>>(d_con, d_listVol, d_list, nActiveBlock);
        check_cuda_error(cudaGetLastError());
        cudaDeviceSynchronize();


        //////////////////////////////////////////////////////////////////
        // 6. update active list
        // read back active volume from the device and add
        // active block to active list on the host memory


        nActiveBlock = 0;
        check_cuda_error( cudaMemcpy(h_listVol, d_listVol, blockNum*sizeof(bool), cudaMemcpyDeviceToHost) );

        for(uint i=0; i<blockNum; i++)
        {
            if(h_listVol[i]) // true : active block (not converged)
            {
                h_listed[i] = true;
                h_list[nActiveBlock++] = i;
            }
            else h_listed[i] = false;
        }
        cudaDeviceSynchronize();

        if (verbose) {
            printf("Iteration : %d\n", nTotalIter);
        }
    }

    if (verbose) {
        printf("Eikonal solver converged after %d iterations\n", nTotalIter);
        printf("Total # of blocks processed : %d\n", nTotalBlockProcessed);
    }

    // delete dynamically allocated host memory
    free(h_list);
    free(h_listed);
    free(h_listVol);
}
