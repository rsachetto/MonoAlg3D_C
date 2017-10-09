//
// Created by sachetto on 05/10/17.
//

//Every model need to implement the functions described in this model file in order to be loaded correctly from the
//edo solver. This models should compile without using any dependency of our codebase

#ifndef MONOALG3D_MODEL_COMMON_H
#define MONOALG3D_MODEL_COMMON_H

#include "../utils/constants.h"

struct cell_model_data {
    int number_of_ode_equations;
    Real initial_v;
    void *extra_data; //TODO: maybe we will need this?
};
#endif //MONOALG3D_MODEL_COMMON_H
