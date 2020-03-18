
#include "utils/file_utils.h"

#ifdef COMPILE_OPENGL
#include "draw/draw.h"
#endif

#include "string/sds.h"
#include "single_file_libraries/stb_ds.h"
#include "vtk_utils/vtk_unstructured_grid.h"
#include "vtk_utils/pvd_utils.h"
#include "draw/draw.h"

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

static void init_draw_config(struct draw_config *draw_config, struct visualization_options *options, bool only_restart) {

    draw_config->grid_info.vtk_grid = NULL;

    draw_config->simulating = true;
    draw_config->exit = false;
    draw_config->restart = false;

    draw_config->paused = true;
    draw_config->advance_or_return = 0;
    draw_config->grid_info.loaded = false;

    if(!only_restart) {
        draw_config->max_v = options->max_v;
        draw_config->min_v = options->min_v;

        if(draw_config->min_v == 0) {
            draw_config->min_v = 0.001f;
        }

        draw_config->dt = options->dt;
        draw_config->draw_type = DRAW_FILE;
        draw_config->grid_info.file_name = NULL;
        draw_config->error_message = NULL;
        draw_config->int_scale = false;
    }

}

static void read_and_render_activation_map(char *input_file) {

    draw_config.grid_info.file_name = NULL;

    omp_set_lock(&draw_config.draw_lock);
    draw_config.grid_info.vtk_grid = new_vtk_unstructured_grid_from_file(input_file);
    draw_config.grid_info.loaded = true;
    draw_config.int_scale = true;

    if(!draw_config.grid_info.vtk_grid) {
        char tmp[4096];
        sprintf(tmp, "%s is not an activation map", input_file);
        draw_config.error_message = strdup(tmp);
        omp_unset_lock(&draw_config.draw_lock);
        return;
    }

    draw_config.grid_info.file_name = input_file;
    draw_config.min_v = draw_config.grid_info.vtk_grid->min_v;
    draw_config.max_v = draw_config.grid_info.vtk_grid->max_v;

    omp_unset_lock(&draw_config.draw_lock);
}

static int read_and_render_files(struct visualization_options *options) {

    bool using_pvd = (options->pvd_file != NULL);

    struct vtk_files *vtk_files;

    const char *input_dir = options->input_folder;
    const char *prefix = options->files_prefix;
    const char *pvd_file = options->pvd_file;

    int current_file = options->start_file;
    int v_step = options->step;

    if(!using_pvd) {
        vtk_files = (struct vtk_files*) malloc(sizeof(struct vtk_files));
        vtk_files->files_list = NULL;
        vtk_files->timesteps = NULL;
        if(input_dir) {
            vtk_files->files_list = list_files_from_dir_sorted(input_dir, prefix);
        }
    }
    else {
        vtk_files = list_files_from_and_timesteps_from_pvd(pvd_file);
    }

    int num_files = arrlen(vtk_files->files_list);
    sds full_path;

    if(!using_pvd) {
        full_path = sdsnew(input_dir);
    }
    else {
        full_path = sdsnew(get_dir_from_path(pvd_file));
    }

    if(!num_files) {
        char tmp[4096];
        if(*full_path) {
            sprintf(tmp, "No simulations file found in %s", full_path);
        }
        else {
            sprintf(tmp, "No path or pvd file provided!");
        }
        fprintf(stderr, "%s\n", tmp);
        draw_config.error_message = strdup(tmp);
        return SIMULATION_FINISHED;

    }
    
    if(current_file > num_files) {
        fprintf(stderr, "[WARN] start_at value (%d) is greater than the number of files (%d). Setting start_at to %d\n", current_file, num_files, num_files);
        current_file = num_files-1;
    }

    int step;
    int step1;
    int step2 = 0;
    int final_step;
    real_cpu dt = 0;

    if(!using_pvd) {
        step1 = get_step_from_filename(vtk_files->files_list[0]);

        if (num_files > 1) {
            step2 = get_step_from_filename(vtk_files->files_list[1]);
            step = step2 - step1;
        } else {
            step = step1;
        }

        final_step = get_step_from_filename(vtk_files->files_list[num_files - 1]);

        dt = draw_config.dt;

        draw_config.step = step;
        if (dt == 0.0) {
            draw_config.final_time = final_step;

        } else {
            draw_config.final_time = final_step * dt;
        }
    }
    else {
        draw_config.final_time = vtk_files->timesteps[num_files-1];
        draw_config.dt = -1; //we don't care about dt here as the timesteps are in the PVD file
    }

    while(true) {

        if(!using_pvd) {
            if (dt == 0) {
                draw_config.time = get_step_from_filename(vtk_files->files_list[current_file]);
            } else {
                draw_config.time = get_step_from_filename(vtk_files->files_list[current_file]) * dt;
            }
        }
        else {
            draw_config.time = vtk_files->timesteps[current_file];
        }

        sdsfree(full_path);

        if(!using_pvd) {
            full_path = sdsnew(input_dir);
        }
        else {
            full_path = sdsnew(get_dir_from_path(pvd_file));
        }

        full_path = sdscat(full_path, "/");

        draw_config.grid_info.file_name = NULL;

        full_path = sdscat(full_path, vtk_files->files_list[current_file]);
        omp_set_lock(&draw_config.draw_lock);
        free_vtk_unstructured_grid(draw_config.grid_info.vtk_grid);
        draw_config.grid_info.vtk_grid = new_vtk_unstructured_grid_from_file(full_path);
        draw_config.grid_info.loaded = true;

        if(!draw_config.grid_info.vtk_grid) {
            char tmp[4096];
            sprintf(tmp, "Decoder not available for file %s", vtk_files->files_list[current_file]);
            fprintf(stderr, "%s\n", tmp);
            draw_config.error_message = strdup(tmp);
            omp_unset_lock(&draw_config.draw_lock);
            arrfree(vtk_files->files_list);
            arrfree(vtk_files->timesteps);
            free(vtk_files);
            sdsfree(full_path);
            return SIMULATION_FINISHED;

        }

        draw_config.grid_info.file_name = full_path;
        omp_unset_lock(&draw_config.draw_lock);

        omp_set_lock(&draw_config.sleep_lock);

        if(draw_config.restart) {
            draw_config.time = 0.0;
            free_vtk_unstructured_grid(draw_config.grid_info.vtk_grid);
            arrfree(vtk_files->files_list);
            arrfree(vtk_files->timesteps);
            free(vtk_files);
            sdsfree(full_path);
            return RESTART_SIMULATION;
        }
        if(draw_config.exit) {
            arrfree(vtk_files->files_list);
            arrfree(vtk_files->timesteps);
            free(vtk_files);
            sdsfree(full_path);
            return END_SIMULATION;
        }

        //TODO: maybe change how we handle advance_return
        if(draw_config.paused) {
            current_file += draw_config.advance_or_return;
            draw_config.advance_or_return = 0;
            if(current_file < 0) current_file++;
            else if(current_file >= num_files) current_file--;

        }
        else {
            current_file += v_step;
            if(current_file >= num_files) {
                current_file -= v_step;
                draw_config.paused = true;
            }

        }
    }
}

int main(int argc, char **argv) {

    struct visualization_options *options = new_visualization_options();

    parse_visualization_options(argc, argv, options);

    omp_set_num_threads(2);

    if(!options->activation_map) {
        if (!options->input_folder) {
            if (!options->pvd_file) {
                fprintf(stderr, "Error!. You have to provide a pvd file or a input folder!\n");
            }
        } else {
            if (options->pvd_file) {
                fprintf(stderr,
                        "Warning!. You provided a pvd file (%s) and an input folder (%s). The input folder will be ignored!\n",
                        options->pvd_file, options->input_folder);
                free(options->input_folder);
                options->input_folder = NULL;

            }
        }
    }

    if(options->save_activation_only) {
        struct vtk_unstructured_grid *vtk_grid = new_vtk_unstructured_grid_from_file(options->activation_map);
        if(!vtk_grid) {
            fprintf(stderr, "Failed to convert %s\n", options->activation_map);
            exit(EXIT_FAILURE);
        }
        sds save_path = sdsnew(options->activation_map);
        save_path = sdscat(save_path, ".vtu");
        save_vtk_unstructured_grid_as_vtu_compressed(vtk_grid, save_path, 6);
        free_vtk_unstructured_grid(vtk_grid);
        sdsfree(save_path);
    }
    else {
        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            {
                omp_init_lock(&draw_config.draw_lock);
                omp_init_lock(&draw_config.sleep_lock);
                init_draw_config(&draw_config, options, false);
                init_and_open_visualization_window(DRAW_FILE);
            }

            #pragma omp section
            {
                if(options->activation_map) {
                    read_and_render_activation_map(options->activation_map);
                }
                else {
                    int result = read_and_render_files(options);

                    while (result == RESTART_SIMULATION || result == SIMULATION_FINISHED) {
                        if (result == RESTART_SIMULATION) {
                            init_draw_config(&draw_config, options, true);
                            result = read_and_render_files(options);
                        }

                        if (draw_config.restart) result = RESTART_SIMULATION;

                        if (result == END_SIMULATION || draw_config.exit) {
                            break;
                        }
                    }
                }

            }
        }
    }

    free_visualization_options(options);
    return EXIT_SUCCESS;
}