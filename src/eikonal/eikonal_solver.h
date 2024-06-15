#ifndef MONOALG3D_EIKONAL_SOLVER_H
#define MONOALG3D_EIKONAL_SOLVER_H

#include "../common_types/common_types.h"
#include "common_def.h"
#include "../alg/cell/cell.h"

struct eikonal_solver {
	bool verbose;
	bool cuda_mem_created;
	size_t width, height, depth, num_active_cells;
	size_t iters_per_block, solver_type;
	real ***speeds;
	real ***answer;
	size_t **seeds;
	bool *** mask;
    struct cell_node **active_cells;
	CUDAMEMSTRUCT memory_struct;
};

struct eikonal_solver * new_eikonal_solver(bool verbose);
void free_eikonal_solver(struct eikonal_solver *solver);
void solve_eikonal(struct eikonal_solver *solver);
void write_alg(struct eikonal_solver *solver, char *filename);

#endif //MONOALG3D_EIKONAL_SOLVER_H
