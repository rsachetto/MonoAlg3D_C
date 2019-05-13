
#include "utils/file_utils.h"

#ifdef COMPILE_OPENGL
#include "draw/draw.h"
#endif

#include "single_file_libraries/stb_ds.h"
#include "vtk_utils/data_utils.h"
#include "vtk_utils/vtk_unstructured_grid.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

static inline int get_step_from_filename(char *filename) {

    char *ia = filename;

    int int_a = 0;

    for(; *ia; ia++) {
        if(isdigit(*ia))
            int_a = int_a*10 + *ia - '0';
    }

    return int_a;

}

void init_draw_config(struct draw_config *draw_config, struct visualization_options *options) {

    draw_config->grid_info.vtk_grid = NULL;

    draw_config->max_v = options->max_v;
    draw_config->min_v = options->min_v;

    if(draw_config->min_v == 0) draw_config->min_v = 0.1f;

    draw_config->simulating = true;

    draw_config->dt = options->dt;

    draw_config->exit = false;
    draw_config->restart = false;

    draw_config->paused = true;
    draw_config->advance_or_return = 0;
    draw_config->draw_type = DRAW_FILE;
    draw_config->grid_info.vtk_grid = NULL;
    draw_config->grid_info.file_name = NULL;
}

static int read_and_render_files(char *input_dir, char* prefix) {

    string_array vtk_file_list  = list_files_from_dir_sorted(input_dir, prefix);

    int num_files = arrlen(vtk_file_list);
    sds full_path = sdsnew(input_dir);

    if(!num_files) {
        fprintf(stderr, "No simulations file found in %s", full_path);
        return SIMULATION_FINISHED;
    }

    int current_file = 0;

    int step1;
    int step2 = 0;

    step1 = get_step_from_filename(vtk_file_list[0]);
    int step;

    if(num_files > 0) {
        step2 = get_step_from_filename(vtk_file_list[1]);
        step = step2 - step1;
    }
    else {
        step = step1;
    }

    int final_step = get_step_from_filename(vtk_file_list[num_files-1]);
    real_cpu dt = draw_config.dt;

    draw_config.step = step;
    draw_config.final_time = final_step*dt;

    while(true) {

        if(draw_config.restart) {
            draw_config.time = 0.0;
            free_vtk_unstructured_grid(draw_config.grid_info.vtk_grid);
            arrfree(vtk_file_list);
            return RESTART_SIMULATION;
        }
        if(draw_config.exit) {
            arrfree(vtk_file_list);
            return END_SIMULATION;
        }

        draw_config.time = get_step_from_filename(vtk_file_list[current_file])*dt;

        sdsfree(full_path);
        full_path = sdsnew(input_dir);
        full_path = sdscat(full_path, "/");

        draw_config.grid_info.file_name = NULL;

        full_path = sdscat(full_path, vtk_file_list[current_file]);
        omp_set_lock(&draw_config.draw_lock);
        free_vtk_unstructured_grid(draw_config.grid_info.vtk_grid);
        draw_config.grid_info.vtk_grid = new_vtk_unstructured_grid_from_vtu_file(full_path);
        draw_config.grid_info.file_name = full_path;
        omp_unset_lock(&draw_config.draw_lock);

        omp_set_lock(&draw_config.sleep_lock);

        if(draw_config.paused) {
            current_file += draw_config.advance_or_return;
            if(current_file < 0) current_file++;
            else if(current_file >= num_files) current_file--;

        }
        else {
            current_file++;
            if(current_file >= num_files) {
                current_file--;
                draw_config.paused = true;
            }

        }
    }

}

int main(int argc, char **argv) {

    struct visualization_options *options = new_visualization_options();

    parse_visualization_options(argc, argv, options);

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
            int result = read_and_render_files(options->input_folder, "V_it_");

            while (result == RESTART_SIMULATION || result == SIMULATION_FINISHED) {
                if(result == RESTART_SIMULATION) {
                    init_draw_config(&draw_config, options);
                    result = read_and_render_files(argv[1], "V_it_");
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
