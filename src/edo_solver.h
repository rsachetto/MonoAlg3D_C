//
// Created by sachetto on 02/10/17.
//

#ifndef MONOALG3D_EDO_SOLVER_H
#define MONOALG3D_EDO_SOLVER_H

struct edo_solver {

    int number_of_equations;
    double initial_v;
    double dt;

    //used for the adaptive time step solver
    double previous_dt;
    double time_new;

    //Use dinamic libraries to load from a .so model file
    //https://www.dwheeler.com/program-library/Program-Library-HOWTO/x172.html
    //https://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html


};

#endif //MONOALG3D_EDO_SOLVER_H
