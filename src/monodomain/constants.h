//
// Created by sachetto on 03/10/17.
//

#ifndef MONOALG3D_CONSTANTS_H
#define MONOALG3D_CONSTANTS_H

#include <math.h>

#define UM2_TO_CM2 0.00000001f
#define BLOCK_SIZE 256

// Precision to be used for the calculations
typedef float real;
typedef double real_cpu;

#define EXPORT_FN

#define ALPHA(beta, cm, dt, dx, dy, dz) ((((beta) * (cm)) / (dt)) * UM2_TO_CM2) * ((dx)*(dy)*(dz))

#define KAPPA(beta, cm, l, h) ((pow(l,4) - pow(h,4)) / (12.0*pow(l,2))) * beta * cm * UM2_TO_CM2 

#endif //MONOALG3D_CONSTANTS_H
