#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "3dparty/sds/sds.h"
#include "3dparty/stb_ds.h"
#include "config/config_parser.h"
#include "gui/gui.h"
#include "utils/file_utils.h"
#include "vtk_utils/pvd_utils.h"
#include "vtk_utils/vtk_unstructured_grid.h"

#define MAX_ERROR_SIZE 4096

static void read_and_render_activation_map(struct gui_shared_info *gui_config, char *input_file, char *error) {

    gui_config->grid_info.file_name = NULL;

    omp_set_lock(&gui_config->draw_lock);
    gui_config->grid_info.vtk_grid = new_vtk_unstructured_grid_from_file(input_file);
    gui_config->grid_info.loaded = true;
    gui_config->int_scale = true;

    if(!gui_config->grid_info.vtk_grid) {
        snprintf(error, MAX_ERROR_SIZE, "%s is not an activation map", input_file);
        if(gui_config->error_message) {
            free(gui_config->error_message);
        }
        gui_config->error_message = strdup(error);
        omp_unset_lock(&gui_config->draw_lock);
        return;
    }

    gui_config->grid_info.file_name = input_file;
    gui_config->min_v = gui_config->grid_info.vtk_grid->min_v;
    gui_config->max_v = gui_config->grid_info.vtk_grid->max_v;

    omp_unset_lock(&gui_config->draw_lock);
}

static void read_visible_cells(struct vtk_unstructured_grid *vtk_grid, sds full_path) {

    sds full_path_cp = sdsnew(full_path);
    full_path_cp = sdscat(full_path_cp, ".vis");
    FILE *vis_file = fopen(full_path_cp, "rw");

    if(vis_file) {
        uint32_t n_cells = vtk_grid->num_cells;
        arrsetlen(vtk_grid->cell_visibility, n_cells);
        fread(vtk_grid->cell_visibility, sizeof(uint8_t), n_cells, vis_file);
        fclose(vis_file);
    }
}

static int read_and_render_files(struct visualization_options *options, struct gui_shared_info *gui_config) {

    char error[MAX_ERROR_SIZE];

    bool using_pvd = false;
    bool single_file = false;
    struct simulation_files *simulation_files = NULL;

    const char *input = options->input;
    const char *prefix = options->files_prefix;
    gui_config->current_file_index = options->start_file;
    int v_step = options->step;

    bool maybe_ensight = false;

    struct path_information input_info;

    get_path_information(input, &input_info);

    if(!input_info.exists) {
        snprintf(error, MAX_ERROR_SIZE,
                 "Invalid path or pvd file provided! Press 'o' to open an directory or 'f' to open a simulation file (pvd, vtu, vtk, acm or alg)!");
        if(gui_config->error_message) {
            free(gui_config->error_message);
        }

        gui_config->error_message = strdup(error);
        return SIMULATION_FINISHED;
    }

    if(input_info.is_dir) {
        simulation_files = (struct simulation_files *)malloc(sizeof(struct simulation_files));
        simulation_files->files_list = NULL;
        simulation_files->timesteps = NULL;
        if(input) {
            string_array ignore_files = NULL;
            arrput(ignore_files, strdup("vis"));
            simulation_files->files_list = list_files_from_dir(input, prefix, NULL, ignore_files, true);

            if(arrlen(simulation_files->files_list) == 0) {
                //Maybe is an ensight folder, lets try again
                simulation_files->files_list = list_files_from_dir(input, "Vm", NULL, ignore_files, true);
                maybe_ensight = true;
            }

            arrfree(ignore_files);
        }
    } else {
        if(strcmp(input_info.file_extension, "pvd") == 0) {
            using_pvd = true;
            simulation_files = list_files_from_and_timesteps_from_pvd(input);
        } else if(strcmp(input_info.file_extension, "acm") == 0) {
            read_and_render_activation_map(gui_config, (char *)input, error);
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

    free_path_information(&input_info);

    uint32_t num_files = arrlen(simulation_files->files_list);

    sds full_path;

    if(!using_pvd) {
        full_path = sdsnew(input);
    } else {
        full_path = sdsnew(get_dir_from_path(input));
    }

    if(!num_files) {
        snprintf(error, MAX_ERROR_SIZE, "No simulations file found in %s", full_path);

        if(gui_config->error_message)
            free(gui_config->error_message);
        gui_config->error_message = strdup(error);

        sdsfree(full_path);
        free(simulation_files);
        return SIMULATION_FINISHED;
    }

    sds geometry_file = NULL;

    if(maybe_ensight) {
        geometry_file = sdscatfmt(sdsempty(), "%s/geometry.geo", input);

        get_path_information(geometry_file,  &input_info);

        if(!input_info.exists) {
            snprintf(error, MAX_ERROR_SIZE, "Geometry file %s not found", geometry_file);

            if(gui_config->error_message)
                free(gui_config->error_message);
            gui_config->error_message = strdup(error);

            sdsfree(full_path);
            free(simulation_files);
            free_path_information(&input_info);

            return SIMULATION_FINISHED;
        }

    }

    if(gui_config->current_file_index > num_files) {
        fprintf(stderr, "[WARN] start_at value (%d) is greater than the number of files (%d). Setting start_at to %d\n", (int)gui_config->current_file_index,
                num_files, num_files);
        gui_config->current_file_index = num_files - 1;
    }

    float dt = 0;
    if(maybe_ensight) {
        //TODO: make this better
        gui_config->dt = 0;
        gui_config->step = 1;
        gui_config->final_file_index = num_files - 1;
        gui_config->final_time = gui_config->final_file_index;
        gui_config->time = -1;

    } else if(!using_pvd) {

        if(single_file) {
            gui_config->dt = -1;
            gui_config->step = 1;
            gui_config->final_file_index = num_files - 1;
            gui_config->final_time = gui_config->final_file_index;
            gui_config->time = -1;
        }
        else {
            int step;
            int step1;
            int final_step;

            step1 = get_step_from_filename(simulation_files->files_list[0]);

            if(num_files > 1) {
                int step2 = 0;
                step2 = get_step_from_filename(simulation_files->files_list[1]);
                step = step2 - step1;
            } else {
                step = step1;
            }

            final_step = get_step_from_filename(simulation_files->files_list[num_files - 1]);

            dt = gui_config->dt;

            gui_config->step = step;

            gui_config->final_file_index = final_step / step;

            if(dt == 0.0) {
                gui_config->final_time = (float)final_step;

            } else {
                gui_config->final_time = (float)final_step * dt;
            }
        }
    } else {
        gui_config->final_time = simulation_files->timesteps[num_files - 1];
        gui_config->dt = -1;
    }

    bool ensigth_grid_loaded = false;

    while(true) {

        omp_set_lock(&gui_config->draw_lock);

        if(maybe_ensight) {
            gui_config->time = (int)gui_config->current_file_index;
        } else if(!using_pvd) {
            if(dt == 0) {
                gui_config->time = (float)get_step_from_filename(simulation_files->files_list[(int)gui_config->current_file_index]);
            } else {
                gui_config->time = (float)get_step_from_filename(simulation_files->files_list[(int)gui_config->current_file_index]) * dt;
            }
        } else {
            gui_config->time = simulation_files->timesteps[(int)gui_config->current_file_index];
        }

        sdsfree(full_path);

        if(!using_pvd) {
            full_path = sdsnew(input);
        } else {
            full_path = sdsnew(get_dir_from_path(input));
        }

        if(!single_file) {
            full_path = sdscat(full_path, "/");
            full_path = sdscat(full_path, simulation_files->files_list[(int)gui_config->current_file_index]);
        }

        if(maybe_ensight) {
            if(!ensigth_grid_loaded) {
                gui_config->grid_info.vtk_grid = new_vtk_unstructured_grid_from_file_with_index(geometry_file, options->value_index);
                ensigth_grid_loaded = true;
            }
            set_vtk_grid_values_from_ensight_file(gui_config->grid_info.vtk_grid, full_path);
        } else {
            free_vtk_unstructured_grid(gui_config->grid_info.vtk_grid);
            gui_config->grid_info.vtk_grid = new_vtk_unstructured_grid_from_file_with_index(full_path, options->value_index);
        }

        if(single_file) {
            gui_config->min_v = gui_config->grid_info.vtk_grid->min_v;
            gui_config->max_v = gui_config->grid_info.vtk_grid->max_v;
        }

        if(!gui_config->grid_info.vtk_grid) {
            snprintf(error, MAX_ERROR_SIZE, "Decoder not available for file %s", simulation_files->files_list[(int)gui_config->current_file_index]);

            if(gui_config->error_message) {
                free(gui_config->error_message);
            }

            gui_config->error_message = strdup(error);
            gui_config->grid_info.loaded = false;
            gui_config->paused = true;
        } else {
            read_visible_cells(gui_config->grid_info.vtk_grid, full_path);
            gui_config->grid_info.file_name = full_path;
            gui_config->grid_info.loaded = true;
        }

        omp_unset_lock(&gui_config->draw_lock);

        // here we wait until the mesh was rendered
        omp_set_lock(&gui_config->sleep_lock);

        if(gui_config->restart) {
            gui_config->time = 0.0f;
            free_vtk_unstructured_grid(gui_config->grid_info.vtk_grid);
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            free(simulation_files);
            sdsfree(full_path);
            return RESTART_SIMULATION;
        }

        if(gui_config->exit) {
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            free(simulation_files);
            sdsfree(full_path);
            return END_SIMULATION;
        }

        if(!gui_config->paused) {
            gui_config->current_file_index += v_step;
            if(gui_config->current_file_index >= num_files) {
                gui_config->current_file_index -= v_step;
                gui_config->paused = true;
            }
        }
    }
}

static void init_gui_config_for_visualization(struct visualization_options *options, struct gui_shared_info *gui_config, bool only_restart) {

    gui_config->grid_info.vtk_grid = NULL;

    gui_config->simulating = true;
    gui_config->exit = false;
    gui_config->restart = false;

    gui_config->paused = true;
    gui_config->grid_info.loaded = false;

    gui_config->ui_scale = options->ui_scale;

    if(!only_restart) {
        gui_config->input = NULL;
        omp_init_lock(&gui_config->draw_lock);
        omp_init_lock(&gui_config->sleep_lock);
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

int main(int argc, char **argv) {

    struct gui_shared_info *gui_config = MALLOC_ONE_TYPE(struct gui_shared_info);

    struct visualization_options *options = new_visualization_options();

    parse_visualization_options(argc, argv, options);

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
                init_gui_config_for_visualization(options, gui_config, false);
                init_and_open_gui_window(gui_config);
            }

            OMP(section) {
                int result = read_and_render_files(options, gui_config);

                while(result == RESTART_SIMULATION || result == SIMULATION_FINISHED) {

                    if(gui_config->input) {
                        options->input = gui_config->input;
                        result = RESTART_SIMULATION;
                    }

                    if(gui_config->restart) {
                        result = RESTART_SIMULATION;
                        gui_config->grid_info.loaded = false;
                    }

                    if(result == RESTART_SIMULATION) {
                        if(options->input) {
                            init_gui_config_for_visualization(options, gui_config, true);
                            result = read_and_render_files(options, gui_config);
                        }
                    } else if(result == END_SIMULATION || gui_config->exit) {
                        break;
                    }
                }
            }
        }
    }

    free_visualization_options(options);
    return EXIT_SUCCESS;
}
