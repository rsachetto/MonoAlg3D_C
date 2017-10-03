//
// Created by sachetto on 02/10/17.
//

#include "ode_solver.h"
#include <stdio.h>
#include <stdlib.h>

void set_ode_initial_conditions(struct ode_solver *solver) {
    if (solver->gpu) {
        if (solver->method == EULER_METHOD) {
//o pitch eh inicializado nesta funcao.
//TODO: @Incomplete
//pitch = setIC_ode_gpu(&sv, originalNumCells);
        } else {
//TODO: @Incomplete
//pitch = setIC_ode_gpu_adapt(&sv, originalNumCells, dt_edo);
        }
    } else {
    //TODO: here we have to call the dinamic loaded funcion
        //setInitalConditionsODEs();
    }
}


const char* get_ODE_method(int met) {

    switch(met) {
        case 0: return "Euler Method";
        case 1: return "Euler Method with adaptive time step (Formula)";
        default: printf("Invalid Method!!\n");  exit(0);
    }
}