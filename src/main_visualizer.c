
#include "utils/file_utils.h"

#include "gui/gui.h"

#include "3dparty/sds/sds.h"
#include "3dparty/stb_ds.h"
#include "config/config_parser.h"
#include "vtk_utils/pvd_utils.h"
#include "vtk_utils/vtk_unstructured_grid.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char error[4096];

static int read_and_render_files(struct visualization_options *options) {

    bool using_pvd = false;
    bool single_file = false;
    struct simulation_files *simulation_files = NULL;

    const char *input = options->input;
    const char *prefix = options->files_prefix;
    int current_file = options->start_file;
    int v_step = options->step;

    struct path_information input_info;

    get_path_information(input, &input_info);

    if(!input_info.exists) {
        sprintf(error, "Invalid path or pvd file provided! Press 'o' to open an directory or 'f' to open a simulation file (pvd, vtu, vtk or acm)!");
        if(gui_get_error_message())
            gui_free_error_message();

        gui_set_error_message(error);
        return SIMULATION_FINISHED;
    }

    if(input_info.is_dir) {
        simulation_files = (struct simulation_files *)malloc(sizeof(struct simulation_files));
        simulation_files->files_list = NULL;
        simulation_files->timesteps = NULL;
        if(input) {
            simulation_files->files_list = list_files_from_dir(input, prefix, NULL, true);
        }
    } else {
        if(strcmp(input_info.file_extension, "pvd") == 0) {
            using_pvd = true;
            simulation_files = list_files_from_and_timesteps_from_pvd(input);
        } else if(strcmp(input_info.file_extension, "acm") == 0) {
            read_and_render_activation_map((char *)input, error);
            return SIMULATION_FINISHED;
        } else if(strcmp(input_info.file_extension, "vtk") == 0 || strcmp(input_info.file_extension, "vtu") == 0 ||
                  strcmp(input_info.file_extension, "txt") == 0 || strcmp(input_info.file_extension, "bin") == 0 ||
                  strcmp(input_info.file_extension, "alg") == 0) {
            simulation_files = (struct simulation_files *)malloc(sizeof(struct simulation_files));
            simulation_files->files_list = NULL;
            simulation_files->timesteps = NULL;
            single_file = true;
            if(input) {
                arrput(simulation_files->files_list, (char *)input);
            }
        }
    }

    int num_files = 0;

    if(simulation_files)
        num_files = arrlen(simulation_files->files_list);

    sds full_path;

    if(!using_pvd) {
        full_path = sdsnew(input);
    } else {
        full_path = sdsnew(get_dir_from_path(input));
    }

    if(!num_files) {
        sprintf(error, "No simulations file found in %s", full_path);

        if(gui_get_error_message())
            gui_free_error_message();
        gui_set_error_message(error);

        sdsfree(full_path);
        free(simulation_files);

        return SIMULATION_FINISHED;
    }

    if(current_file > num_files) {
        fprintf(stderr, "[WARN] start_at value (%d) is greater than the number of files (%d). Setting start_at to %d\n", current_file, num_files, num_files);
        current_file = num_files - 1;
    }

    int step;
    int step1;
    int step2 = 0;
    int final_step;
    real_cpu dt = 0;

    if(!using_pvd) {
        step1 = get_step_from_filename(simulation_files->files_list[0]);

        if(num_files > 1) {
            step2 = get_step_from_filename(simulation_files->files_list[1]);
            step = step2 - step1;
        } else {
            step = step1;
        }

        final_step = get_step_from_filename(simulation_files->files_list[num_files - 1]);

        dt = gui_get_dt();

        //gui_config.step = step;
        gui_set_step(step);
        if(dt == 0.0) {
            gui_set_final_time(final_step);

        } else {
            gui_set_final_time(final_step * dt);
        }
    } else {
        gui_set_final_time(simulation_files->timesteps[num_files - 1]);
        gui_set_dt(-1);
    }

    while(true) {

        if(!using_pvd) {
            if(dt == 0) {
                gui_set_time(get_step_from_filename(simulation_files->files_list[current_file]));
            } else {
                gui_set_time(get_step_from_filename(simulation_files->files_list[current_file]) * dt);
            }
        } else {
            gui_set_time(simulation_files->timesteps[current_file]);
        }

        sdsfree(full_path);

        if(!using_pvd) {
            full_path = sdsnew(input);
        } else {
            if(using_pvd)
                full_path = sdsnew(get_dir_from_path(input));
        }

        gui_set_filename(NULL);

        if(!single_file) {
            full_path = sdscat(full_path, "/");
            full_path = sdscat(full_path, simulation_files->files_list[current_file]);
        }

        gui_lock_draw_lock();

        free_vtk_unstructured_grid(gui_get_vtk_grid());
        gui_set_vtk_grid(new_vtk_unstructured_grid_from_file(full_path));

        gui_set_grid_loaded(true);

        if(!gui_get_vtk_grid()) {
            sprintf(error, "Decoder not available for file %s", simulation_files->files_list[current_file]);
            if(gui_get_error_message())
                gui_free_error_message();

            gui_set_error_message(error);

            gui_unlock_draw_lock();

            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            free(simulation_files);
            sdsfree(full_path);
            return SIMULATION_FINISHED;
        }

        gui_set_filename(full_path);
        gui_unlock_draw_lock();

        gui_unlock_sleep_lock();

        if(gui_get_restart()) {
           gui_set_time(0.0);
            free_vtk_unstructured_grid(gui_get_vtk_grid());
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            free(simulation_files);
            sdsfree(full_path);
            return RESTART_SIMULATION;
        }

        if(gui_get_exit()) {
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            free(simulation_files);
            sdsfree(full_path);
            return END_SIMULATION;
        }

        // TODO: maybe change how we handle advance_return
        if(gui_get_paused()) {
            current_file += gui_get_advance_or_return();
            gui_set_advance_or_return(0);
            if(current_file < 0)
                current_file++;
            else if(current_file >= num_files)
                current_file--;

        } else {
            current_file += v_step;
            if(current_file >= num_files) {
                current_file -= v_step;
                gui_set_paused(true);
            }
        }
    }
}

int main(int argc, char **argv) {

    struct visualization_options *options = new_visualization_options();

    parse_visualization_options(argc, argv, options);

    omp_set_num_threads(2);

    if(!options->input) {
        fprintf(stderr, "[ERR] You have to provide a pvd file, a input folder or an activation map!\n");
    }

    if(options->save_activation_only) {
        struct vtk_unstructured_grid *vtk_grid = new_vtk_unstructured_grid_from_file(options->input);
        if(!vtk_grid) {
            fprintf(stderr, "Failed to convert %s\n", options->input);
            exit(EXIT_FAILURE);
        }
        sds save_path = sdsnew(options->input);
        save_path = sdscat(save_path, ".vtu");
        save_vtk_unstructured_grid_as_vtu_compressed(vtk_grid, save_path, 6);
        free_vtk_unstructured_grid(vtk_grid);
        sdsfree(save_path);
    } else {
        OMP(parallel sections num_threads(2)) {
            OMP(section) {
                init_gui_config_for_visualization(options, false);
                init_and_open_gui_window();
            }

            OMP(section) {
                int result = read_and_render_files(options);

                while(result == RESTART_SIMULATION || result == SIMULATION_FINISHED) {

                    if(gui_get_input()) {
                        options->input = gui_get_input();
                        result = RESTART_SIMULATION;
                    }

                    if(gui_get_restart()) {
                        result = RESTART_SIMULATION;
                        gui_set_grid_loaded(false);
                        //gui_config.grid_info.loaded = false;
                    }

                    if(result == RESTART_SIMULATION) {
                        init_gui_config_for_visualization(options, false);
                        result = read_and_render_files(options);
                    }

                    else if(result == END_SIMULATION || gui_get_exit()) {
                        break;
                    }
                }
            }
        }
    }

    free_visualization_options(options);
    return EXIT_SUCCESS;
}
