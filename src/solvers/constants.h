//
// Created by sachetto on 03/10/17.
//

#ifndef MONOALG3D_CONSTANTS_H
#define MONOALG3D_CONSTANTS_H

#include <math.h>

#define UM2_TO_CM2 0.00000001
#define BLOCK_SIZE 256


// Precision to be used for the calculations
typedef float Real;

inline double ALPHA(double beta, double cm, double dt, double h) {
    return  ( ((beta * cm)/dt ) * UM2_TO_CM2 ) * pow(h,3.0);

}

#endif //MONOALG3D_CONSTANTS_H
