//
// Created by sachetto on 15/03/22.
//

#ifndef MONOALG3D_ECG_CONFIG_H
#define MONOALG3D_ECG_CONFIG_H

#include "../alg/grid/grid.h"
#include "../common_types/common_types.h"
#include "../monodomain/constants.h"
#include "../monodomain/monodomain_solver.h"
#include "config_common.h"

#define CALC_ECG(name) void name(struct time_info *time_info, struct config *config, struct grid *the_grid)
typedef CALC_ECG(calc_ecg_fn);

#define INIT_CALC_ECG(name) void name(struct config *config, struct ode_solver *the_ode_solver, struct grid *the_grid)
typedef INIT_CALC_ECG(init_calc_ecg_fn);

#define END_CALC_ECG(name) void name(struct config *config)
typedef END_CALC_ECG(end_calc_ecg_fn);

#define CALL_INIT_CALC_ECG(config, ode_solver, grid)                                                                                                           \
    do {                                                                                                                                                       \
        if((config) && (config)->init_function) {                                                                                                              \
            ((init_calc_ecg_fn *)(config)->init_function)((config), (ode_solver), (grid));                                                                           \
        }                                                                                                                                                      \
    } while(0)

#define CALL_END_CALC_ECG(config)                                                                                                                              \
    do {                                                                                                                                                       \
        if((config) && (config)->end_function) {                                                                                                               \
            ((end_calc_ecg_fn *)(config)->end_function)((config));                                                                                               \
        }                                                                                                                                                      \
    } while(0)

#define print_calc_ecg_config_values(s) LOG_COMMON_CONFIG("[ecg]", s)

#endif // MONOALG3D_ECG_CONFIG_H
