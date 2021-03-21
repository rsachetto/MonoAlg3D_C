//
// Created by sachetto on 19/03/2021.
//

#include "helper_functions.h"
#include <stdlib.h>
#include "../config_helpers/config_helpers.h"


struct extra_data_for_fibrosis* set_common_schemia_data(struct config *config, uint32_t num_cells) {

    struct extra_data_for_fibrosis *extra_data = MALLOC_ONE_TYPE(struct extra_data_for_fibrosis);

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

    extra_data->atpi = atpi;
    extra_data->Ko = Ko;
    extra_data->Ki = Ki;
    extra_data->GNa_multiplicator = GNa_multiplicator;
    extra_data->GCaL_multiplicator = GCaL_multiplicator;
    extra_data->INaCa_multiplicator = INaCa_multiplicator;
    extra_data->Vm_modifier = Vm_modifier;
    extra_data->fibrosis = MALLOC_ARRAY_OF_TYPE(real, num_cells);

    return extra_data;
}