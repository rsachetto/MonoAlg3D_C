//
// Created by sachetto on 11/11/17.
//

#ifndef MONOALG3D_DRAW_H
#define MONOALG3D_DRAW_H

#include "../alg/grid/grid.h"
#include "../alg/cell/cell.h"
#include "../monodomain/ode_solver.h"
#include "../vtk_utils/vtk_unstructured_grid.h"

#include <omp.h>

#define DRAW_SIMULATION 0
#define DRAW_FILE 1

#define MIN_VERTICAL_TICKS 4
#define MAX_VERTICAL_TICKS 20

#define MIN_HORIZONTAL_TICKS 4
#define MAX_HORIZONTAL_TICKS 20

struct action_potential {
    real_cpu v;
    real_cpu t;
};

typedef struct action_potential * action_potential_array;

struct draw_config {

    real_cpu max_v;
    real_cpu min_v;
    bool simulating;
    bool paused;
    bool exit;
    bool adaptive;
    bool restart;
    real_cpu time;
    real_cpu final_time;
    real_cpu dt;
    int step;

    char *config_name;

    long solver_time;
    long ode_total_time;
    long cg_total_time;
    long total_mat_time;
    long total_ref_time;
    long total_deref_time;
    long total_write_time;
    long total_config_time;
    long total_cg_it;

    int advance_or_return;
    int draw_type;
    //If we are compiling this file, openmp is available.
    omp_lock_t draw_lock;
    omp_lock_t sleep_lock;

    struct grid_info {
        union {
            struct vtk_unstructured_grid *vtk_grid;
            struct grid *grid_to_draw;
        };
        char *file_name;
        bool loaded;
    } grid_info;

    char *error_message;

};

struct draw_config draw_config;

void init_and_open_visualization_window();

#endif //MONOALG3D_DRAW_H
