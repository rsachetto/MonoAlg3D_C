//
// Created by sachetto on 19/03/2021.
//

#include "helper_functions.h"
#include <stdlib.h>
#include "../config_helpers/config_helpers.h"


real* set_common_schemia_data(struct config *config, uint32_t num_cells, int num_par, size_t *extra_data_size) {

    *extra_data_size = sizeof(real)*(num_cells + num_par);

    real *extra_data = (real*)malloc(*extra_data_size);

    real atpi = 6.8;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, atpi, config, "atpi");

    real Ko = 5.4;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ko, config, "Ko");

    real Ki = 138.3;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ki, config, "Ki");

    real GNa_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GNa_multiplicator, config, "GNa_multiplicator");

    real GCaL_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GCaL_multiplicator, config, "GCaL_multiplicator");

    real INaCa_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaCa_multiplicator, config, "INaCa_multiplicator");

    real Vm_modifier = 0.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Vm_modifier, config, "Vm_modifier");

    extra_data[0] = atpi;
    extra_data[1] = Ko;
    extra_data[2] = Ki;
    extra_data[3] = Vm_modifier;
    extra_data[4] = GNa_multiplicator;
    extra_data[5] = GCaL_multiplicator;
    extra_data[6] = INaCa_multiplicator;

    return extra_data;
}