//
// CUDA implementation of FIM (Fast Iterative Method) for Eikonal equations
//
// Copyright (c) Won-Ki Jeong (wkjeong@unist.ac.kr)
//
// 2016. 2. 4
//

#include "cuda_fim_kernel.h"

__device__ real get_time_eikonal(real a, real b, real c, real s)
{
	real ret, tmp;

	// a > b > c
	if(a < b) { tmp = a; a = b; b = tmp; }
	if(b < c) { tmp = b; b = c; c = tmp; }
	if(a < b) { tmp = a; a = b; b = tmp; }

	ret = INF;

	if(c < INF)
	{
		ret = c + s;
		
		if(ret > b) 
		{	
			tmp = ((b+c) + sqrtf(2.0f*s*s-(b-c)*(b-c)))*0.5f;
		
			if(tmp > b) ret = tmp; 

			if(ret > a)	{			
				tmp = (a+b+c)/3.0f + sqrtf(2.0f*(a*(b-a)+b*(c-b)+c*(a-c))+3.0f*s*s)/3.0f; 

				if(tmp > a) ret = tmp;
			}
		}
	}

	return ret;
}

__global__ void run_solver(real* spd, bool* mask, const real *sol_in, real *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, int nIter, uint nActiveBlock)
{
	uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

	if(list_idx < nActiveBlock)
	{
		// retrieve actual block index from the active list
		uint block_idx = list[list_idx];

		real F;
		bool isValid;
		uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH;
		uint base_addr = block_idx*blocksize;

		uint xgridlength = xdim/BLOCK_LENGTH;
		uint ygridlength = ydim/BLOCK_LENGTH;
		uint zgridlength = zdim/BLOCK_LENGTH;

		// compute block index
		uint bx = block_idx%xgridlength;
		uint tmpIdx = (block_idx - bx)/xgridlength;
		uint by = tmpIdx%ygridlength;
		uint bz = (tmpIdx-by)/ygridlength;

		uint tx = threadIdx.x;
		uint ty = threadIdx.y;
		uint tz = threadIdx.z;
		uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

		__shared__ real _sol[BLOCK_LENGTH+2][BLOCK_LENGTH+2][BLOCK_LENGTH+2];

		// copy global to shared memory
		dim3 idx(tx+1,ty+1,tz+1);

		SOL(idx.x,idx.y,idx.z) = sol_in[base_addr + tIdx];
		F = spd[base_addr + tIdx];
		if(F > 0) F = 1.0/F; // F = 1/f
		isValid = mask[base_addr + tIdx];

		uint new_base_addr, new_tIdx;

		// 1-neighborhood values
		if(tx == 0) 
		{
			if(bx == 0) // end of the grid
			{	
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx + BLOCK_LENGTH-1;
				new_base_addr = (block_idx - 1)*blocksize;	
			}

			SOL(tx,idx.y,idx.z) = sol_in[new_base_addr + new_tIdx];	
		}

		if(tx == BLOCK_LENGTH-1)
		{
			if(bx == xgridlength-1) // end of the grid
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx - (BLOCK_LENGTH-1);
				new_base_addr = (block_idx + 1)*blocksize;	
			}
			SOL(tx+2,idx.y,idx.z) = sol_in[new_base_addr + new_tIdx];	
		}

		if(ty == 0)
		{
			if(by == 0)
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH;
				new_base_addr = (block_idx - xgridlength)*blocksize;
			}

			SOL(idx.x,ty,idx.z) = sol_in[new_base_addr + new_tIdx];
		}

		if(ty == BLOCK_LENGTH-1)
		{
			if(by == ygridlength-1) 
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH;
				new_base_addr = (block_idx + xgridlength)*blocksize;
			}

			SOL(idx.x,ty+2,idx.z) = sol_in[new_base_addr + new_tIdx];
		}

		if(tz == 0)
		{
			if(bz == 0)
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
				new_base_addr = (block_idx - xgridlength*ygridlength)*blocksize;
			}

			SOL(idx.x,idx.y,tz) = sol_in[new_base_addr + new_tIdx];
		}

		if(tz == BLOCK_LENGTH-1)
		{
			if(bz == zgridlength-1) 
			{
				new_tIdx = tIdx;
				new_base_addr = base_addr;
			}
			else
			{
				new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
				new_base_addr = (block_idx + xgridlength*ygridlength)*blocksize;
			}

			SOL(idx.x,idx.y,tz+2) = sol_in[new_base_addr + new_tIdx];
		}

		__syncthreads();

		real a,b,c,oldT,newT;

		for(int iter=0; iter<nIter; iter++)	
		{
			//
			// compute new value
			//
			oldT = newT = SOL(idx.x,idx.y,idx.z);

			if(isValid)
			{
				a = min(SOL(tx,idx.y,idx.z),SOL(tx+2,idx.y,idx.z));
				b = min(SOL(idx.x,ty,idx.z),SOL(idx.x,ty+2,idx.z));
				c = min(SOL(idx.x,idx.y,tz),SOL(idx.x,idx.y,tz+2));

				real tmp = (real) get_time_eikonal(a, b, c, F);

				newT = min(tmp,oldT);
			}
			__syncthreads();	

			if(isValid) SOL(idx.x,idx.y,idx.z) = newT;

			__syncthreads(); // this may not required		
		}


		real residue = oldT - newT;

		// write back to global memory
		con[base_addr + tIdx] = (residue < EPS) ? true : false;
		sol_out[base_addr + tIdx] = newT;		
	}
}

__global__ void run_reduction(bool *con, bool *listVol, uint *list, uint nActiveBlock)
{
	uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

	if(list_idx < nActiveBlock)
	{
		uint block_idx = list[list_idx];

		__shared__ bool conv[BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH];

		uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH/2;
		uint base_addr = block_idx*blocksize*2;
		uint tx = threadIdx.x;
		uint ty = threadIdx.y;
		uint tz = threadIdx.z;
		uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

		conv[tIdx] = con[base_addr + tIdx];
		conv[tIdx + blocksize] = con[base_addr + tIdx + blocksize];

		__syncthreads();

		for(uint i=blocksize; i>0; i/=2)
		{
			if(tIdx < i)
			{
				bool b1, b2;
				b1 = conv[tIdx];
				b2 = conv[tIdx+i];
				conv[tIdx] = (b1 && b2) ? true : false ;
			}
			__syncthreads();
		}

		if(tIdx == 0) 
		{		
			listVol[block_idx] = !conv[0]; // active list is negation of tile convergence (active = not converged)
		}
	}
}

__global__ void run_check_neighbor(real* spd, bool* mask, const real *sol_in, real *sol_out, 
																	 bool *con, uint* list, int xdim, int ydim, int zdim,
																	 uint nActiveBlock, uint nTotalBlock)
{

	uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

	if(list_idx < nTotalBlock)
	{
		real F;
		bool isValid;
		__shared__ real _sol[BLOCK_LENGTH+2][BLOCK_LENGTH+2][BLOCK_LENGTH+2];

		uint block_idx = list[list_idx];
		uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH;
		uint base_addr = block_idx*blocksize;

		uint tx = threadIdx.x;
		uint ty = threadIdx.y;
		uint tz = threadIdx.z;
		uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

		if(list_idx < nActiveBlock) // copy value
		{
			sol_out[base_addr + tIdx] = sol_in[base_addr + tIdx];
		} 
		else
		{
			uint xgridlength = xdim/BLOCK_LENGTH;
			uint ygridlength = ydim/BLOCK_LENGTH;
			uint zgridlength = zdim/BLOCK_LENGTH;

			// compute block index
			uint bx = block_idx%xgridlength;
			uint tmpIdx = (block_idx - bx)/xgridlength;
			uint by = tmpIdx%ygridlength;
			uint bz = (tmpIdx-by)/ygridlength;

			// copy global to shared memory
			dim3 idx(tx+1,ty+1,tz+1);
			_sol[idx.x][idx.y][idx.z] = sol_in[base_addr + tIdx];
			F = spd[base_addr + tIdx];
			if(F > 0) F = 1.0/F;
			isValid = mask[base_addr + tIdx];

			uint new_base_addr, new_tIdx;

			// 1-neighborhood values
			if(tx == 0) 
			{
				if(bx == 0) // end of the grid
				{	
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx + BLOCK_LENGTH-1;
					new_base_addr = (block_idx - 1)*blocksize;	
				}
				_sol[tx][idx.y][idx.z] = sol_in[new_base_addr + new_tIdx];	
			}

			if(tx == BLOCK_LENGTH-1)
			{
				if(bx == xgridlength-1) // end of the grid
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx - (BLOCK_LENGTH-1);
					new_base_addr = (block_idx + 1)*blocksize;	
				}
				_sol[tx+2][idx.y][idx.z] = sol_in[new_base_addr + new_tIdx];	
			}

			if(ty == 0)
			{
				if(by == 0)
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH;
					new_base_addr = (block_idx - xgridlength)*blocksize;
				}
				_sol[idx.x][ty][idx.z] = sol_in[new_base_addr + new_tIdx];
			}

			if(ty == BLOCK_LENGTH-1)
			{
				if(by == ygridlength-1) 
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH;
					new_base_addr = (block_idx + xgridlength)*blocksize;
				}
				_sol[idx.x][ty+2][idx.z] = sol_in[new_base_addr + new_tIdx];
			}

			if(tz == 0)
			{
				if(bz == 0)
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
					new_base_addr = (block_idx - xgridlength*ygridlength)*blocksize;
				}
				_sol[idx.x][idx.y][tz] = sol_in[new_base_addr + new_tIdx];
			}

			if(tz == BLOCK_LENGTH-1)
			{
				if(bz == zgridlength-1) 
				{
					new_tIdx = tIdx;
					new_base_addr = base_addr;
				}
				else
				{
					new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
					new_base_addr = (block_idx + xgridlength*ygridlength)*blocksize;
				}
				_sol[idx.x][idx.y][tz+2] = sol_in[new_base_addr + new_tIdx];
			}

			__syncthreads();


			real a,b,c,oldT,newT;

			//
			// compute new value
			//
			oldT = newT = _sol[idx.x][idx.y][idx.z];

			if(isValid)
			{
				a = min(_sol[tx][idx.y][idx.z],_sol[tx+2][idx.y][idx.z]);
				b = min(_sol[idx.x][ty][idx.z],_sol[idx.x][ty+2][idx.z]);
				c = min(_sol[idx.x][idx.y][tz],_sol[idx.x][idx.y][tz+2]);

				real tmp = (real) get_time_eikonal(a, b, c, F);
				newT = min(tmp,oldT);

				sol_out[base_addr + tIdx] = newT;
			}
			// write back to global memory
			real residue = oldT - newT;
			con[base_addr + tIdx] = (residue < EPS) ? true : false;	
		}
	}
}


