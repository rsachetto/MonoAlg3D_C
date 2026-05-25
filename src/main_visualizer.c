#include <ctype.h>
#include <float.h>
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
#include <unistd.h>

#define MAX_ERROR_SIZE 4096

static void read_and_render_activation_map(struct gui_shared_info *gui_config, char *input_file, char *error) {

    gui_config->grid_info.file_name = NULL;

    omp_set_nest_lock(&gui_config->draw_lock);
    gui_config->grid_info.vtk_grid = new_vtk_unstructured_grid_from_file(input_file, true);
    gui_config->grid_info.loaded = true;
    gui_config->int_scale = true;

    if(!gui_config->grid_info.vtk_grid) {
        snprintf(error, MAX_ERROR_SIZE, "%s is not an activation map", input_file);
        if(gui_config->message) {
            free(gui_config->message);
        }
        gui_config->message = strdup(error);
        omp_unset_nest_lock(&gui_config->draw_lock);
        return;
    }

    gui_config->grid_info.file_name = input_file;
    gui_config->min_v = gui_config->grid_info.vtk_grid->min_v;
    gui_config->max_v = gui_config->grid_info.vtk_grid->max_v;

    omp_unset_nest_lock(&gui_config->draw_lock);
}

static void calc_vm_bounds(struct gui_shared_info *gui_config, const struct simulation_files *simulation_files, const sds geometry_file, const bool ensight) {

    char path[2048];
    gui_config->grid_info.loaded = false;

    char error[MAX_ERROR_SIZE];
    snprintf(error, MAX_ERROR_SIZE, "Calculating Vm bounds!");

    if(gui_config->message) {
        free(gui_config->message);
    }

    gui_config->message = strdup(error);
    static struct vtk_unstructured_grid *tmp_grid = NULL;
    static bool en_tmp_loaded = false;

    gui_config->max_v = -10000.0;
    gui_config->min_v = 10000.0;

    uint32_t num_files = arrlen(simulation_files->files_list);
    gui_config->file_size = num_files;

    for(uint32_t i = 0; i < num_files; i++) {
        gui_config->progress = i;
        const char *current_file_name = simulation_files->files_list[i];
        sprintf(path, "%s/%s", simulation_files->base_dir, current_file_name);

        if(ensight) {
            if(!en_tmp_loaded) {
                tmp_grid = new_vtk_unstructured_grid_from_file(geometry_file, true);
                en_tmp_loaded = true;
            }

            set_vtk_grid_values_from_ensight_file(tmp_grid, path);
        } else {
            tmp_grid = new_vtk_unstructured_grid_from_file(path, true);
        }

        if(tmp_grid->max_v > gui_config->max_v) {
            gui_config->max_v = tmp_grid->max_v;
        }

        if(tmp_grid->min_v < gui_config->min_v) {
            gui_config->min_v = tmp_grid->min_v;
        }

        if(!ensight) {
            free_vtk_unstructured_grid(tmp_grid);
        }
    }

    if(ensight) {
        free_vtk_unstructured_grid(tmp_grid);
        tmp_grid = NULL;
        en_tmp_loaded = false;
    }

    gui_config->grid_info.loaded = true;
    gui_config->calc_bounds = false;
}

static int read_and_render_files(struct visualization_options *options, struct gui_shared_info *gui_config) {

    char error[MAX_ERROR_SIZE];

    bool using_pvd = false;
    bool single_file = false;
    struct simulation_files *simulation_files = NULL;

    const char *input = options->input;
    const char *prefix = options->files_prefix;
    gui_config->current_file_index = (float)options->start_file;
    int v_step = options->step;

    bool ensight = false;

    struct path_information input_info;

    get_path_information(input, &input_info);
    bool esca_file = false;
    bool geo_file = false;

    if(!input_info.exists) {
        snprintf(error, MAX_ERROR_SIZE,
                 "Invalid path or pvd file provided! Press 'o' to open an directory or 'f' to open a simulation file (pvd, vtu, vtk, acm or alg)!");
        if(gui_config->message) {
            free(gui_config->message);
        }

        gui_config->message = strdup(error);
        gui_config->input = NULL;
        options->input = NULL;
        return SIMULATION_FINISHED;
    }

    if(input_info.is_dir) {
        simulation_files = (struct simulation_files *)malloc(sizeof(struct simulation_files));
        simulation_files->files_list = NULL;
        simulation_files->timesteps = NULL;

        string_array ignore_files = NULL;
        arrput(ignore_files, strdup("vis"));
        simulation_files->files_list = list_files_from_dir(input, prefix, NULL, ignore_files, true);

        if(arrlen(simulation_files->files_list) == 0) {
            // Maybe is an ensight folder, lets try again
            simulation_files->files_list = list_files_from_dir(input, "Vm", NULL, ignore_files, true);
            ensight = true;
        }

        arrfree(ignore_files);
    } else {
        if(FILE_HAS_EXTENSION(input_info, "pvd")) {
            using_pvd = true;
            simulation_files = list_files_from_and_timesteps_from_pvd(input);
        } else if(FILE_HAS_EXTENSION(input_info, "acm")) {
            read_and_render_activation_map(gui_config, (char *)input, error);
            return SIMULATION_FINISHED;
        } else if(FILE_HAS_EXTENSION(input_info, "vtk") || FILE_HAS_EXTENSION(input_info, "vtu") || FILE_HAS_EXTENSION(input_info, "txt") ||
                  FILE_HAS_EXTENSION(input_info, "bin") || FILE_HAS_EXTENSION(input_info, "alg") || FILE_HAS_EXTENSION(input_info, "geo") ||
                  FILE_HAS_EXTENSION_PREFIX(input_info, "Esca")) {
            simulation_files = (struct simulation_files *)malloc(sizeof(struct simulation_files));
            simulation_files->files_list = NULL;
            simulation_files->timesteps = NULL;
            single_file = true;

            if(FILE_HAS_EXTENSION(input_info, "geo") || FILE_HAS_EXTENSION_PREFIX(input_info, "Esca")) {
                ensight = true;
                if(FILE_HAS_EXTENSION(input_info, "geo"))
                    geo_file = true;
                else
                    esca_file = true;
            }

            arrput(simulation_files->files_list, (char *)input);
        }
    }

    if(!using_pvd) {
        simulation_files->base_dir = sdsnew(input);
    } else {
        simulation_files->base_dir = sdsnew(get_dir_from_path(input));
    }

    uint32_t num_files = arrlen(simulation_files->files_list);

    if(!num_files) {
        snprintf(error, MAX_ERROR_SIZE, "No simulations file found in %s", simulation_files->base_dir);

        if(gui_config->message)
            free(gui_config->message);
        gui_config->message = strdup(error);

        sdsfree(simulation_files->base_dir);
        free(simulation_files);
        return SIMULATION_FINISHED;
    }

    sds geometry_file = NULL;

    if(ensight) {
        if(!single_file || esca_file) {
            if(esca_file) {
                geometry_file = sdscatfmt(sdsempty(), "%s/geometry.geo", input_info.dir_name);
            } else {
                geometry_file = sdscatfmt(sdsempty(), "%s/geometry.geo", input);
            }
            get_path_information(geometry_file, &input_info);
        } else {
            geometry_file = sdsnew(input);
        }

        if(!input_info.exists) {
            snprintf(error, MAX_ERROR_SIZE, "Geometry file %s not found", geometry_file);

            if(gui_config->message)
                free(gui_config->message);
            gui_config->message = strdup(error);

            sdsfree(simulation_files->base_dir);
            free(simulation_files);
            free_path_information(&input_info);

            return SIMULATION_FINISHED;
        }
    }

    if(gui_config->current_file_index > (float)num_files) {
        fprintf(stderr, "[WARN] start_at value (%d) is greater than the number of files (%u). Setting start_at to %u\n", (int)gui_config->current_file_index,
                num_files, num_files);
        gui_config->current_file_index = (float)(num_files - 1);
    }

    float dt = 0;
    if(!using_pvd) {

        if(single_file) {

            gui_config->dt = -1;
            gui_config->step = 1;
            gui_config->final_file_index = (int)num_files - 1;
            gui_config->final_time = (float)gui_config->final_file_index;

        } else if(ensight) {

            struct path_information case_info;
            sds case_file_path = sdscatfmt(sdsempty(), "%s/simulation_result.case", input);
            get_path_information(case_file_path, &case_info);

            if(!case_info.exists) {
                gui_config->dt = 1;
            } else {
                simulation_files->timesteps = read_timesteps_from_case_file(case_file_path);
            }
            sdsfree(case_file_path);
            free_path_information(&case_info);

            if(arrlen(simulation_files->files_list) != arrlen(simulation_files->timesteps)) {
                printf("[WARN] Inconsistent number of files and file steps\n");
            }

            gui_config->dt = -1;
            gui_config->final_file_index = (int)num_files - 1;
            gui_config->final_time = simulation_files->timesteps[num_files - 1];

        } else {

            int step;
            int step1;
            int final_step;

            step1 = get_step_from_filename(simulation_files->files_list[0]);

            if(num_files > 1) {
                int step2 = get_step_from_filename(simulation_files->files_list[1]);
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
        gui_config->final_file_index = (int)num_files - 1;
        gui_config->final_time = simulation_files->timesteps[num_files - 1];
        gui_config->dt = -1;
    }

    bool ensigth_grid_loaded = false;
    bool ensigth_vis_loaded = false;

    free_path_information(&input_info);

    gui_config->simulation_files = simulation_files;

    char full_path[2048];

    if(ensight || single_file) {
        gui_config->enable_slice = true;
    }

    while(true) {

        if(!single_file && gui_config->calc_bounds) {
            calc_vm_bounds(gui_config, simulation_files, geometry_file, ensight);
        }

        int current_file_index = (int)gui_config->current_file_index;
        char *current_file_name = simulation_files->files_list[current_file_index];

        if(single_file) {
            gui_config->time = gui_config->current_file_index;
            sprintf(full_path, "%s", simulation_files->base_dir);
        } else {
            if(using_pvd || ensight) {
                gui_config->time = simulation_files->timesteps[current_file_index];
            } else {
                gui_config->time = (float)get_step_from_filename(current_file_name);
                if(dt != 0) {
                    gui_config->time = gui_config->time * dt;
                }
            }
            sprintf(full_path, "%s/%s", simulation_files->base_dir, current_file_name);
        }

        omp_set_nest_lock(&gui_config->draw_lock);

        if(ensight) {
            if(!ensigth_grid_loaded) {
                gui_config->grid_info.vtk_grid = new_vtk_unstructured_grid_from_file(geometry_file, single_file);
                ensigth_grid_loaded = true;
            }

            if(!geo_file) {
                set_vtk_grid_values_from_ensight_file(gui_config->grid_info.vtk_grid, full_path);
            } else {
                gui_config->grid_info.vtk_grid->min_v = 0.0001f;
                gui_config->grid_info.vtk_grid->max_v = 0.0002f;
            }

        } else {
            free_vtk_unstructured_grid(gui_config->grid_info.vtk_grid);
            gui_config->grid_info.vtk_grid =
                new_vtk_unstructured_grid_from_file_with_progress(full_path, single_file, &gui_config->progress, &gui_config->file_size);
        }

        if(!gui_config->grid_info.vtk_grid) {
            snprintf(error, MAX_ERROR_SIZE, "Decoder not available for file %s", current_file_name);

            if(gui_config->message) {
                free(gui_config->message);
            }

            gui_config->message = strdup(error);
            gui_config->grid_info.loaded = false;
            gui_config->paused = true;
        } else {
            if(ensight) {
                if(!ensigth_vis_loaded) {
                    read_or_calc_visible_cells(&gui_config->grid_info.vtk_grid, geometry_file);
                    ensigth_vis_loaded = true;
                }
            } else {
                read_or_calc_visible_cells(&gui_config->grid_info.vtk_grid, full_path);
            }
            // TODO: for ensigth, maybe we should put the data name here.
            gui_config->grid_info.file_name = full_path;
            gui_config->grid_info.loaded = true;
        }

        if(single_file) {
            gui_config->max_v = gui_config->grid_info.vtk_grid->max_v;
            gui_config->min_v = gui_config->grid_info.vtk_grid->min_v;
        }

        omp_unset_nest_lock(&gui_config->draw_lock);

        // here we wait until the mesh was rendered
        omp_set_nest_lock(&gui_config->sleep_lock);

        if(gui_config->restart) {
            gui_config->time = 0.0f;
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            sdsfree(simulation_files->base_dir);
            free(simulation_files);
            return RESTART_SIMULATION;
        }

        if(gui_config->exit) {
            arrfree(simulation_files->files_list);
            arrfree(simulation_files->timesteps);
            sdsfree(simulation_files->base_dir);
            free(simulation_files);
            return END_SIMULATION;
        }

        if(!gui_config->paused) {
            gui_config->current_file_index += (float)v_step;
            if(gui_config->current_file_index >= (float)num_files) {
                gui_config->current_file_index -= (float)v_step;
                gui_config->paused = true;
            }
        }
    }
}

static void init_gui_config_for_visualization(const struct visualization_options *options, struct gui_shared_info *gui_config, bool only_restart) {

    // TODO: set this from command line
    gui_config->adaptive = false;

    gui_config->grid_info.vtk_grid = NULL;

    gui_config->simulating = true;
    gui_config->exit = false;
    gui_config->restart = false;

    gui_config->paused = true;
    gui_config->grid_info.loaded = false;

    gui_config->enable_slice = false;

    gui_config->ui_scale = options->ui_scale;

    if(!only_restart) {
        gui_config->input = NULL;
        omp_init_nest_lock(&gui_config->draw_lock);
        omp_init_nest_lock(&gui_config->sleep_lock);
        gui_config->max_v = options->max_v;
        gui_config->min_v = options->min_v;

        if(gui_config->min_v == 0) {
            gui_config->min_v = 0.001f;
        }

        gui_config->dt = options->dt;
        gui_config->draw_type = DRAW_FILE;
        gui_config->grid_info.file_name = NULL;
        gui_config->message = NULL;
        gui_config->int_scale = false;
    }
}

int main(int argc, char **argv) {

    struct gui_shared_info *gui_config = CALLOC_ONE_TYPE(struct gui_shared_info);

    struct visualization_options *options = new_visualization_options();

    parse_visualization_options(argc, argv, options);

    if(options->save_activation_only) {
        struct vtk_unstructured_grid *vtk_grid = new_vtk_unstructured_grid_from_file(options->input, false);
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

                    // HACK: this should not be needed, we have to find a way to avoid this hack.
                    // If we take this out, the open mesh option does not work properly
                    usleep(10);

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
                    } else if(gui_config->exit) {
                        break;
                    }
                }
            }
        }
    }

    free_visualization_options(options);
    return EXIT_SUCCESS;
}
