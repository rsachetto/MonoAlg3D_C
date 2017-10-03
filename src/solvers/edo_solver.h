//
// Created by sachetto on 02/10/17.
//

#ifndef MONOALG3D_EDO_SOLVER_H
#define MONOALG3D_EDO_SOLVER_H

#include <stdbool.h>

#define EULER_METHOD 0
#define EULER_METHOD_ADPT 1

struct edo_solver {

    int number_of_equations;
    double initial_v;
    double min_dt;
    int method;

    //used for the adaptive time step solver
    double previous_dt;
    double time_new;

    bool gpu;
    int gpu_id;

    float *sv;

    //Use dinamic libraries to load from a .so model file
    //https://www.dwheeler.com/program-library/Program-Library-HOWTO/x172.html
    //https://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html


};

void set_edo_initial_conditions(struct edo_solver *solver);

#endif //MONOALG3D_EDO_SOLVER_H
