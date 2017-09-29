//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_SOLVER_H
#define MONOALG3D_SOLVER_H

#include <unitypes.h>

struct solver {
    int numberOfIterations;
    double pError;
    double roundStimSize;
    uint8_t parallel;
    uint8_t gpu;
    int method;

    // Time used for solving wave equation.
    double delta_t;
};

#endif //MONOALG3D_SOLVER_H
