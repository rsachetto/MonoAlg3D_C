
#include "utils/file_utils.h"

#ifdef COMPILE_OPENGL
#include "draw/draw.h"
#endif

#include "yxml.h"
#include "single_file_libraries/stb_ds.h"
#include "vtk_utils/data_utils.h"
#include "vtk_utils/vtk_unstructured_grid.h"

#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <yxml.h>

#include <ctype.h>

void init_draw_config(struct draw_config *draw_config, struct user_options *options) {

    draw_config->vtk_grid = NULL;

    draw_config->max_v = options->max_v;
    draw_config->min_v = options->min_v;

    if(draw_config->min_v == 0) draw_config->min_v = 0.1f;

    draw_config->simulating = false;
    draw_config->time = 0.0;

    draw_config->adaptive = options->adaptive;
    draw_config->final_time = options->final_time;
    draw_config->dt = options->dt_pde;

    draw_config->exit = false;
    draw_config->restart = false;
}

struct vtk_unstructured_grid *vtk_grid;

static int read_files(char *input_dir, char* prefix) {

    string_array vtk_file_list  = list_files_from_dir_ordered(input_dir, prefix);

    int num_files = arrlen(vtk_file_list);

    for(int f = 0; f < num_files; f++) {

        sds full_path = sdsnew(input_dir);
        full_path = sdscat(full_path, "/");

        if(draw_config.restart) {
            draw_config.time = 0.0;
            return RESTART_SIMULATION;
        }
        if(draw_config.exit) return END_SIMULATION;

        if(!draw_config.paused) {
            omp_set_lock(&draw_config.draw_lock);
            full_path = sdscat(full_path, vtk_file_list[f]);
            vtk_grid = new_vtk_unstructured_grid_from_vtu_file(full_path);
            draw_config.vtk_grid = vtk_grid;
            omp_unset_lock(&draw_config.draw_lock);
            if(f == 0) draw_config.paused = true;
        }
        else {
            //Is this a good usage of mutexes to mimic sleep-wakeup???
            omp_set_lock(&draw_config.sleep_lock);
            continue;
        }

        //TODO: how to free this memory??
        //free_vtk_unstructured_grid(vtk_grid);
        sdsfree(full_path);
    }

    return SIMULATION_FINISHED;


}

int main(int argc, char **argv) {

    struct user_options *options = new_user_options();

    //TODO: parse command line for this parameters
    options->max_v = 40.0;
    options->min_v = -86;

    options->adaptive = true;
    options->final_time = 0.0;
    options->dt_pde = 0.02;

    draw_config.paused = false;
    //struct draw_config->

#pragma omp parallel sections num_threads(2)
    {
#pragma omp section
        {
            omp_init_lock(&draw_config.draw_lock);
            omp_init_lock(&draw_config.sleep_lock);
            init_draw_config(&draw_config, options);
            init_and_open_visualization_window(DRAW_FILE);
        }

#pragma omp section
        {
            int result = read_files(argv[1], "V_it_");

            while (result == RESTART_SIMULATION || result == SIMULATION_FINISHED) {
                if(result == RESTART_SIMULATION) {
                    init_draw_config(&draw_config, options);
                    result = read_files(argv[1], "V_it_");
                }

                if(draw_config.restart) result = RESTART_SIMULATION;

                if(draw_config.exit)  {
                    break;
                }
            }

        }
    }
    return EXIT_SUCCESS;
}
