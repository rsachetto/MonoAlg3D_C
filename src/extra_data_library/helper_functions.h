//
// Created by sachetto on 19/03/2021.
//

#ifndef MONOALG3D_C_EXTRA_DATA_HELPER_FUNCTIONS_H
#define MONOALG3D_C_EXTRA_DATA_HELPER_FUNCTIONS_H

#include <unistd.h>
#include "../common_types/common_types.h"
#include "../config/config_common.h"

#define SET_EXTRA_DATA_SIZE(value) *extra_data_size = (value)

struct extra_data_for_fibrosis {
    real atpi;
    real Ko;
    real Ki;
    real GNa_multiplicator;
    real GCaL_multiplicator;
    real INaCa_multiplicator;
    real Vm_modifier;
    real *fibrosis;
};

struct extra_data_for_fibrosis * set_common_schemia_data(struct config *config, uint32_t num_cells);

#endif // MONOALG3D_C_EXTRA_DATA_HELPER_FUNCTIONS_H
