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

#ifdef _MSC_VER
#define EXPORT_FN __declspec(dllexport)
#else
#define EXPORT_FN 
#endif

#define ALPHA(beta, cm, dt, dx, dy, dz) ((((beta) * (cm)) / (dt)) * UM2_TO_CM2) * ((dx)*(dy)*(dz))

// ***********************************************************************************************************
// Berg and Pedro code ...
#define ALPHA_CM(beta, cm, dt, dx, dy, dz) (((1.0) / (dt))) * ((dx)*(dy)*(dz))
// ***********************************************************************************************************

#endif //MONOALG3D_CONSTANTS_H
