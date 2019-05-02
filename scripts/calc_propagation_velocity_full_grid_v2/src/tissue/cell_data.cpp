#include "cell_data.h"

struct cell_data* new_cell_data (const uint32_t total_num_cells, const uint32_t total_num_steps)
{
    struct cell_data *result = (struct cell_data*)malloc(sizeof(struct cell_data)*total_num_steps);

    result->num_steps = total_num_steps;

    for (uint32_t i = 0; i < total_num_steps; i++)
        result[i].vms = (double*)malloc(sizeof(double)*total_num_cells);

    return result;
}

void free_cell_data (struct cell_data *data)
{
    if (data)
    {
        free(data->vms);
        free(data);
    }
}