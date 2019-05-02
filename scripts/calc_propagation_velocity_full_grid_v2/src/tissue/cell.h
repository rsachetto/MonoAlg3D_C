#ifndef CALC_PROPAGATION_VELOCITY_CELL_H
#define CALC_PROPAGATION_VELOCITY_CELL_H

#include <iostream>
#include <string>

#include <cstdio>
#include <cstdlib>

struct cell
{
    uint32_t id;                // Tissue index

    double x, y, z;             // Coordinates
    double at;                  // Activation time
    double cv;                  // Conduction velocity

    bool in_boundary_1;         // East
    bool in_boundary_2;         // South
    bool in_boundary_3;         // West
    bool in_boundary_4;         // North
};

void set_new_cell (struct cell *the_cell, const uint32_t id,\
                        const double x, const double y, const double z);
void set_boundaries (struct cell *the_cell, const uint32_t i, const uint32_t j,\
            const uint32_t nx, const uint32_t ny);


#endif //MONOALG3D_UTILS_H_H