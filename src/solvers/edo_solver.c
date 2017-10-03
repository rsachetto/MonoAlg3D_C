//
// Created by sachetto on 02/10/17.
//

#include "edo_solver.h"

void set_edo_initial_conditions(struct edo_solver *solver) {
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