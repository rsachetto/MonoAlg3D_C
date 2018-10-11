//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_RESTORE_STATE_CONFIG_H
#define MONOALG3D_RESTORE_STATE_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../monodomain/ode_solver.h"
#include "../alg/grid/grid.h"

//Forward declaration
struct restore_state_config;

#define RESTORE_STATE(name) EXPORT_FN void name(char* output_dir,                    \
                                                struct restore_state_config *config, \
                                                struct grid *the_grid,               \
                                                struct ode_solver *the_ode_solver)
typedef RESTORE_STATE(restore_state_fn);


struct restore_state_config {
    struct config_common config_data;    
    char *in_dir_name;
    bool in_dir_name_was_set;

    restore_state_fn *restore_state;     
};

struct restore_state_config* new_save_restore_config();
void init_restore_state_functions(struct restore_state_config *config);
void free_restore_state_config(struct restore_state_config* s);
void print_restore_state_config_values(struct restore_state_config* s);


#endif //MONOALG3D_SAVE_CONFIG_H
