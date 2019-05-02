#include "cell.h"

void set_new_cell (struct cell *the_cell, const uint32_t id,\
                        const double x, const double y, const double z)
{
    the_cell->id = id;
    
    the_cell->x = x;
    the_cell->y = y;
    the_cell->z = z;

    the_cell->at = 0.0;
    the_cell->cv = 0.0;

    the_cell->in_boundary_1 = false;
    the_cell->in_boundary_2 = false;
    the_cell->in_boundary_3 = false;
    the_cell->in_boundary_4 = false;
}

void set_boundaries (struct cell *the_cell, const uint32_t i, const uint32_t j,\
            const uint32_t nx, const uint32_t ny)
{
    // Inside east boundary
    if (i == 0)
        the_cell->in_boundary_1 = true;
    // Inside south boundary
    if (j == (ny-1))
        the_cell->in_boundary_2 = true;
    // Inside west boundary
    if (i == (nx-1))
        the_cell->in_boundary_3 = true;
    // Inside north boundary
    if (j == 0)
        the_cell->in_boundary_4 = true;
}