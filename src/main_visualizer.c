
#include "utils/file_utils.h"

#include "gui/gui.h"

#include "3dparty/sds/sds.h"
#include "3dparty/stb_ds.h"
#include "vtk_utils/vtk_unstructured_grid.h"
#include "vtk_utils/pvd_utils.h"
#include "config/config_parser.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

static char error[4096];

static inline int get_step_from_filename(char *filename) {
    char *ia = filename;

    int int_a = 0;

    for(; *ia; ia++) {
        if(isdigit(*ia))
            int_a = int_a*10 + *ia - '0';
    }

    return int_a;
}

static void init_gui_config(struct gui_config *gui_config, struct visualization_options *options, bool only_restart) {

    gui_config->grid_info.vtk_grid = NULL;

    gui_config->simulating = true;
    gui_config->exit = false;
    gui_config->restart = false;

    gui_config->paused = true;
    gui_config->advance_or_return = 0;
    gui_config->grid_info.loaded = false;
    gui_config->input = NULL;

    if(!only_restart) {
        gui_config->max_v = options->max_v;
        gui_config->min_v = options->min_v;

        if(gui_config->min_v == 0) {
            gui_config->min_v = 0.001f;
        }

        gui_config->dt = options->dt;
        gui_config->draw_type = DRAW_FILE;
        gui_config->grid_info.file_name = NULL;
        gui_config->error_message = NULL;
        gui_config->int_scale = false;
    }

}

static void read_and_render_activation_map(char *input_file) {

    gui_config.grid_info.file_name = NULL;

    omp_set_lock(&gui_config.draw_lock);
    gui_config.grid_info.vtk_grid = new_vtk_unstructured_grid_from_file(input_file);
    gui_config.grid_info.loaded = true;
    gui_config.int_scale = true;

    if(!gui_config.grid_info.vtk_grid) {
        sprintf(error, "%s is not an activation map", input_file);
        if(gui_config.error_message) free(gui_config.error_message);
        gui_config.error_message = strdup(error);
        omp_unset_lock(&gui_config.draw_lock);
        return;
    }

    gui_config.grid_info.file_name = input_file;
    gui_config.min_v = gui_config.grid_info.vtk_grid->min_v;
    gui_config.max_v = gui_config.grid_info.vtk_grid->max_v;

    omp_unset_lock(&gui_config.draw_lock);
}

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
        if(gui_config.error_message) free(gui_config.error_message);
        gui_config.error_message = strdup(error);
        return SIMULATION_FINISHED;
    }

    if(input_info.is_dir) {
        simulation_files = (struct simulation_files*) malloc(sizeof(struct simulation_files));
        simulation_files->files_list = NULL;
        simulation_files->timesteps = NULL;
        if(input) {
            simulation_files->files_list = list_files_from_dir_sorted(input, prefix);
        }
    }
    else {
        if(strcmp(input_info.file_extension, "pvd") == 0) {
            using_pvd = true;
            simulation_files = list_files_from_and_timesteps_from_pvd(input);
        } else if(strcmp(input_info.file_extension, "acm") == 0) {
            read_and_render_activation_map((char *)input);
            return SIMULATION_FINISHED;
        } else if(strcmp(input_info.file_extension, "vtk") == 0 || strcmp(input_info.file_extension, "vtu") == 0 ||  strcmp(input_info.file_extension, "txt") == 0
            || strcmp(input_info.file_extension, "bin") == 0 || strcmp(input_info.file_extension, "alg") == 0) {
            simulation_files = (struct simulation_files *)malloc(sizeof(struct simulation_files));
            simulation_files->files_list = NULL;
            simulation_files->timesteps = NULL;
            single_file = true;
            if(input) {
                arrput(simulation_files->files_list, (char*)input);
            }
        }
    }

    int num_files = 0;

    if(simulation_files) num_files = arrlen(simulation_files->files_list);

    sds full_path;

    if(!using_pvd) {
        full_path = sdsnew(input);
    }
    else {
        full_path = sdsnew(get_dir_from_path(input));
    }

    if(!num_files) {
        sprintf(error, "No simulations file found in %s", full_path);

        if(gui_config.error_message) free(gui_config.error_message);
        gui_config.error_message = strdup(error);

        sdsfree(full_path);
        free(simulation_files);

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
        step1 = get_step_from_filename(simulation_files->files_list[0]);

        if (num_files > 1) {
            step2 = get_step_from_filename(simulation_files->files_list[1]);
            step = step2 - step1;
        } else {
            step = step1;
        }

        final_step = get_step_from_filename(simulation_files->files_list[num_files - 1]);

        dt = gui_config.dt;

        gui_config.step = step;
        if (dt == 0.0) {
            gui_config.final_time = final_step;

        } else {
            gui_config.final_time = final_step * dt;
        }
    }
    else {
        gui_config.final_time = simulation_files->timesteps[num_files-1];
        gui_config.dt = -1; //we don't care about dt here as the time steps are in the PVD file
    }

    while(true) {

        if(!using_pvd) {
            if (dt == 0) {
                gui_config.time = get_step_from_filename(simulation_files->files_list[current_file]);
            } else {
                gui_config.time = get_step_from_filename(simulation_files->files_list[current_file]) * dt;
            }
        }
        else {
            gui_config.time = simulation_files->timesteps[current_file];
        }

        sdsfree(full_path);

        if(!using_pvd) {
            full_path = sdsnew(input);
        }
        else {
            if(using_pvd)
                full_path = sdsnew(get_dir_from_path(input));
        }

        gui_config.grid_info.file_name = NULL;

        if(!single_file) {
            full_path = sdscat(full_path, "/");
            full_path = sdscat(full_path, simulation_files->files_list[current_file]);
        }

        omp_set_lock(&gui_config.draw_lock);
        free_vtk_unstructured_grid(gui_config.grid_info.vtk_grid);
        gui_config.grid_info.vtk_grid = new_vtk_unstructured_grid_from_file(full_path);
        gui_config.grid_info.loaded = true;

        if(!gui_config.grid_info.vtk_grid) {
            sprintf(error, "Decoder not available for file %s", simulation_files->files_list[current_file]);
            if(gui_config.error_message) free(gui_config.error_message);
            gui_config.error_message = strdup(error);
            omp_unset_lock(&gui_config.draw_lock);
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            free(simulation_files);
            sdsfree(full_path);
            return SIMULATION_FINISHED;
        }

        gui_config.grid_info.file_name = full_path;
        omp_unset_lock(&gui_config.draw_lock);

        omp_set_lock(&gui_config.sleep_lock);

        if(gui_config.restart) {
            gui_config.time = 0.0;
            free_vtk_unstructured_grid(gui_config.grid_info.vtk_grid);
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            free(simulation_files);
            sdsfree(full_path);
            return RESTART_SIMULATION;
        }

        if(gui_config.exit) {
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            free(simulation_files);
            sdsfree(full_path);
            return END_SIMULATION;
        }

        //TODO: maybe change how we handle advance_return
        if(gui_config.paused) {
            current_file += gui_config.advance_or_return;
            gui_config.advance_or_return = 0;
            if(current_file < 0) current_file++;
            else if(current_file >= num_files) current_file--;

        }
        else {
            current_file += v_step;
            if(current_file >= num_files) {
                current_file -= v_step;
                gui_config.paused = true;
            }

        }
    }
}

int main(int argc, char **argv) {

    struct visualization_options *options = new_visualization_options();

    parse_visualization_options(argc, argv, options);

    omp_set_num_threads(2);

    if (!options->input) {
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
    }
    else {
        OMP(parallel sections num_threads(2))
        {
            OMP(section)
            {
                omp_init_lock(&gui_config.draw_lock);
                omp_init_lock(&gui_config.sleep_lock);
                init_gui_config(&gui_config, options, false);
                init_and_open_gui_window();
            }

            OMP(section)
            {
                int result = read_and_render_files(options);

                while (result == RESTART_SIMULATION || result == SIMULATION_FINISHED) {

                    if(gui_config.input) {
                        options->input = gui_config.input;
                        result = RESTART_SIMULATION;
                    }

                    if (gui_config.restart) {
                        result = RESTART_SIMULATION;
                        gui_config.grid_info.loaded = false;
                    }

                    if (result == RESTART_SIMULATION) {
                        init_gui_config(&gui_config, options, true);
                        result = read_and_render_files(options);
                    }

                    else if (result == END_SIMULATION || gui_config.exit) {
                        break;
                    }
                }

            }
        }
    }

    free_visualization_options(options);
    return EXIT_SUCCESS;
}
