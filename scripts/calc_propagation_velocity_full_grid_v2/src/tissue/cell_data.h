#ifndef CALC_PROPAGATION_VELOCITY_CELL_DATA_H
#define CALC_PROPAGATION_VELOCITY_CELL_DATA_H

#include <iostream>
#include <string>

#include <cstdio>
#include <cstdlib>

struct cell_data
{
    uint32_t num_steps;         // Number of timesteps

    double *vms;
};

struct cell_data* new_cell_data (const uint32_t total_num_cells, const uint32_t total_num_steps);
void free_cell_data (struct cell_data *data);


#endif //MONOALG3D_UTILS_H_H