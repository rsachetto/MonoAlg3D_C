//
// Created by sachetto on 11/11/17.
//

#include <bits/stdint-uintn.h>
#include <float.h>
#include <string.h>
#include <sys/types.h>

// This need to be before RAYGUI_IMPLEMENTATION
#include "gui.h"
#include "gui_colors.h"
#include "gui_draw.h"
#include "gui_mesh_helpers.h"
#include "gui_window_helpers.h"

#include "../utils/file_utils.h"

#include "../3dparty/raylib/src/external/glad.h"
#include "../3dparty/raylib/src/rcamera.h"
#include "../3dparty/raylib/src/rlgl.h"
#include "../3dparty/stb_ds.h"
#include "../3dparty/tinyfiledialogs/tinyfiledialogs.h"

#define RAYGUI_IMPLEMENTATION
#include "../3dparty/raylib/src/extras/raygui.h"

#define GUI_TEXTBOX_EXTENDED_IMPLEMENTATION
#include "../3dparty/raylib/src/extras/gui_textbox_extended.h"

#define RLIGHTS_IMPLEMENTATION
#include "../3dparty/raylib/src/extras/rlights.h"

#include "raylib_ext.h"

static const char *help_box_strings[] = {"  Mouse Wheel to Zoom in-out",
                                         "  Mouse Wheel Pressed to Pan",
                                         "  Alt + Mouse Wheel Pressed to Rotate",
                                         "  Alt + Ctrl + Mouse Wheel Pressed for Smooth Zoom",
                                         "  Ctrl + F to search a cell based on it's center",
                                         "  Hold Ctrl and move the mouse over a cell to see it's position",
                                         "  G to only draw the grid lines",
                                         "  L to enable or disable the grid lines",
                                         "  R to restart simulation (only works when paused)",
                                         "  Alt + R to restart simulation and the box positions",
                                         "  X to show/hide AP visualization",
                                         "  Q to show/hide scale",
                                         "  C to show/hide everything except grid",
                                         "  F to open a simulation file",
                                         "  O to open a simulation directory",
                                         "  F12 to take a screenshot",
                                         "  CTRL + F12 to start/stop recording the screen",
                                         "  . or , to change color scales",
                                         "  Right arrow to advance one dt when paused",
                                         "  Hold up arrow to advance time when paused",
                                         "  Double click on a volume to show the AP",
                                         "  from 1 to 0 to show different file column (or PGUP/PGDOWN)",
                                         "  S to enter mesh slice mode",
                                         "  Space to start or pause simulation"};

#define BOX_MARGIN 25.0f

static const char *slice_help_box_strings[] = {" Press backspace to reset and exit slice mode", " Press enter to accept the sliced mesh",
                                               " Move the slicing plane with arrow keys", " Rotate the slicing plane with ALT + arrow keys"};

static const char *edit_help_box_strings[] = {" Press E to exit edit mode", " Press S to save the current file (overwrite)",
                                              " Right click to copy a cell property", "  Left click to paste the copied property to another cell"};

static struct gui_state *new_gui_state_with_font_sizes(float font_size_small, float font_size_big, float ui_scale) {

    struct gui_state *gui_state = CALLOC_ONE_TYPE(struct gui_state);
    gui_state->current_selected_volumes = NULL;

    gui_state->ui_scale = ui_scale;

    set_camera_params(&(gui_state->camera), true);

    gui_state->font = LoadFont("res/Roboto.ttf");

    gui_state->font_size_small = font_size_small * ui_scale;
    gui_state->font_size_big = font_size_big * ui_scale;

    gui_state->font_spacing_big = gui_state->font_size_big / (float)gui_state->font.baseSize;
    gui_state->font_spacing_small = gui_state->font_size_small / (float)gui_state->font.baseSize;

    gui_state->handle_keyboard_input = true;

    gui_state->scale.window.show = true;

    gui_state->help_box.window.show = false;
    gui_state->slice_help_box.window.show = false;
    gui_state->edit_help_box.window.show = false;

    gui_state->ctrl_pressed = false;

    gui_state->current_data_index = -1;

    gui_state->end_info_box.window.show = true;
    gui_state->mesh_info_box.window.show = true;
    gui_state->ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->current_mouse_over_volume.position_draw = (Vector3){-1, -1, -1};
    gui_state->found_volume.position_mesh = (Vector3){-1, -1, -1};
    gui_state->mouse_pos = (Vector2){0, 0};
    gui_state->current_scale = 0;
    gui_state->voxel_alpha = 255;
    gui_state->scale_alpha = 255;
    gui_state->mouse_timer = -1;
    gui_state->search_window.move = false;

    gui_state->ap_graph_config = MALLOC_ONE_TYPE(struct ap_graph_config);
    gui_state->ap_graph_config->selected_ap_point = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd1 = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd2 = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_aps = NULL;
    gui_state->ap_graph_config->graph.drag = false;
    gui_state->ap_graph_config->graph.move = false;
    gui_state->ap_graph_config->graph.show = true;

    gui_state->current_window_width = GetScreenWidth();
    gui_state->current_window_height = GetScreenHeight();

    gui_state->ap_graph_config->graph.bounds.height = 300.0f * ui_scale;
    gui_state->ap_graph_config->graph.bounds.width = 690.0f * ui_scale;

    gui_state->ap_graph_config->graph.bounds.x = 10;
    gui_state->ap_graph_config->graph.bounds.y = (float)gui_state->current_window_height - gui_state->ap_graph_config->graph.bounds.height - 90;

    hmdefault(gui_state->ap_graph_config->selected_aps, NULL);

    gui_state->search_window.bounds.width = 220;
    gui_state->search_window.bounds.height = 100;

    gui_state->scale.window.bounds.x = (float)gui_state->current_window_width - 30.0f * ui_scale;
    gui_state->scale.window.bounds.y = (float)gui_state->current_window_height / 1.5f;

    gui_state->scale.window.bounds.width = 20;
    gui_state->scale.window.bounds.height = 0;
    gui_state->scale.calc_bounds = true;

    gui_state->controls_window.show = true;

    gui_state->current_mode = VISUALIZING;
    gui_state->slicing_mesh = false;
    gui_state->recalculating_visibility = false;

    gui_state->visibility_recalculated = false;
    gui_state->old_cell_visibility = NULL;
    gui_state->exclude_from_mesh = NULL;

    gui_state->plane_roll = 0.0f;
    gui_state->plane_pitch = 0.0f;
    gui_state->plane_tx = 0.0f;
    gui_state->plane_ty = 0.0f;
    gui_state->plane_tz = 0.0f;

    gui_state->mesh_scale_factor = 1.0f;
    gui_state->mesh_offset = (Vector3){0, 0, 0};

    gui_state->show_coordinates = true;
    gui_state->double_clicked = false;

    gui_state->get_cell_property = false;
    gui_state->paste_cell_property = false;

    const u_int8_t NUM_BUTTONS = 7;
    gui_state->controls_window.bounds.width = NUM_BUTTONS * 32.0 + (NUM_BUTTONS + 1) * 4 + 96;
    gui_state->controls_window.bounds.height = 38.0f + WINDOW_STATUSBAR_HEIGHT;

    return gui_state;
}

static inline void configure_end_info_box_strings(struct gui_shared_info *gui_config, char ***info_string) {

    char tmp[TMP_SIZE];

    int index = 0;

    snprintf(tmp, TMP_SIZE, "Resolution Time: %ld s", gui_config->solver_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "ODE Total Time: %ld s", gui_config->ode_total_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "CG Total Time: %ld s", gui_config->cg_total_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Mat time: %ld s", gui_config->total_mat_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Refine time: %ld s", gui_config->total_ref_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Unrefine time: %ld s", gui_config->total_deref_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Write time: %ld s", gui_config->total_write_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "CG Total Iterations: %ld", gui_config->total_cg_it);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Final Time: %.3lf ms", gui_config->time);
    (*(info_string))[index] = strdup(tmp);
}

static inline bool configure_mesh_info_box_strings(struct gui_state *gui_state, struct gui_shared_info *gui_config, char ***info_string, int draw_type,
                                                   struct mesh_info *mesh_info) {

    if(!gui_config->grid_info.alg_grid && !gui_config->grid_info.vtk_grid)
        return false;

    char tmp[TMP_SIZE];

    uint32_t n_active = 0;

    if(draw_type == DRAW_SIMULATION) {
        n_active = gui_config->grid_info.alg_grid->num_active_cells;
    } else {
        n_active = gui_config->grid_info.vtk_grid->num_cells;
    }

    snprintf(tmp, TMP_SIZE, "Num. of Volumes: %u", n_active);
    arrput((*(info_string)), strdup(tmp));

    snprintf(tmp, TMP_SIZE, "Max X: %f", mesh_info->max_size.x);
    arrput((*(info_string)), strdup(tmp));

    snprintf(tmp, TMP_SIZE, "Max Y: %f", mesh_info->max_size.y);
    arrput((*(info_string)), strdup(tmp));

    snprintf(tmp, TMP_SIZE, "Max Z: %f", mesh_info->max_size.z);
    arrput((*(info_string)), strdup(tmp));

    snprintf(tmp, TMP_SIZE, "Min X: %f", mesh_info->min_size.x);
    arrput((*(info_string)), strdup(tmp));

    snprintf(tmp, TMP_SIZE, "Min Y: %f", mesh_info->min_size.y);
    arrput((*(info_string)), strdup(tmp));

    snprintf(tmp, TMP_SIZE, "Min Z: %f", mesh_info->min_size.z);
    arrput((*(info_string)), strdup(tmp));

    if(!gui_config->adaptive) {
        if(draw_type == DRAW_FILE) {
            snprintf(tmp, TMP_SIZE, "Average DX: %f", gui_config->grid_info.vtk_grid->average_discretization.x);
            arrput((*(info_string)), strdup(tmp));

            snprintf(tmp, TMP_SIZE, "Average DY: %f", gui_config->grid_info.vtk_grid->average_discretization.y);
            arrput((*(info_string)), strdup(tmp));

            snprintf(tmp, TMP_SIZE, "Average DZ: %f", gui_config->grid_info.vtk_grid->average_discretization.z);
            arrput((*(info_string)), strdup(tmp));
        }
    }

    if(draw_type == DRAW_SIMULATION) {

        snprintf(tmp, TMP_SIZE, "--");
        arrput((*(info_string)), strdup(tmp));

        if(gui_config->paused) {
            snprintf(tmp, TMP_SIZE, "Simulation paused: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
        } else if(gui_config->simulating) {
            snprintf(tmp, TMP_SIZE, "Simulation running: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
        } else {
            snprintf(tmp, TMP_SIZE, "Simulation finished: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
            gui_config->paused = true;
        }

        arrput((*(info_string)), strdup(tmp));

    } else {
        if(gui_state->max_data_index > 0) {
            snprintf(tmp, TMP_SIZE, "Current column idx: %d of %d", gui_state->current_data_index + 2, gui_state->max_data_index + 1);
        } else {
            snprintf(tmp, TMP_SIZE, "--");
        }

        arrput((*(info_string)), strdup(tmp));

        if(gui_config->final_file_index == 0) {
            snprintf(tmp, TMP_SIZE, "Visualizing single file.                      ");
        } else {
            if(gui_config->paused) {
                if(gui_config->dt == 0) {
                    snprintf(tmp, TMP_SIZE, "Visualization paused: %d of %d steps", (int)gui_config->time, (int)gui_config->final_time);
                } else {
                    snprintf(tmp, TMP_SIZE, "Visualization paused: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
                }
            } else if(gui_config->simulating) {
                if(gui_config->dt == 0) {
                    snprintf(tmp, TMP_SIZE, "Visualization running: %d of %d steps", (int)gui_config->time, (int)gui_config->final_time);
                } else {
                    snprintf(tmp, TMP_SIZE, "Visualization running: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
                }

            } else {
                if(gui_config->dt == 0) {
                    snprintf(tmp, TMP_SIZE, "Visualization finished: %d of %d steps", (int)gui_config->time, (int)gui_config->final_time);
                } else {
                    snprintf(tmp, TMP_SIZE, "Visualization finished: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
                }
            }
        }

        arrput((*(info_string)), strdup(tmp));
    }

    Vector2 wider_text = MeasureTextEx(gui_state->font, tmp, gui_state->font_size_small, gui_state->font_spacing_small);

    float box_w = wider_text.x * 1.08f;

    gui_state->mesh_info_box.window.bounds.width = box_w;

    return true;
}

static void handle_keyboard_input(struct gui_shared_info *gui_config, struct gui_state *gui_state) {

    int kp = GetKeyPressed();

    if(gui_config->draw_type == DRAW_FILE && gui_config->final_file_index == 0) {
        // Visualizing single file
        if(kp == KEY_E) {
            if(gui_state->current_mode == VISUALIZING) {
                gui_state->current_mode = EDITING;
            } else if(gui_state->current_mode == EDITING) {
                gui_state->current_mode = VISUALIZING;
            }
        }
    }

    if(gui_config->paused) {

        if(IsKeyUp(KEY_RIGHT_CONTROL) || IsKeyUp(KEY_LEFT_CONTROL)) {
            gui_state->ctrl_pressed = false;
        }

        if(IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown(KEY_LEFT_CONTROL)) {

            gui_state->ray_mouse_over = GetMouseRay(GetMousePosition(), gui_state->camera);

            gui_state->ctrl_pressed = true;

            if(kp == KEY_S) {
                char const *filter[1] = {"*.vtu"};

                const char *save_path = tinyfd_saveFileDialog("Save VTK file", gui_config->input, 1, filter, "vtu files");

                if(save_path) {
                    save_vtk_unstructured_grid_as_vtu_compressed(gui_config->grid_info.vtk_grid, save_path, 6);
                    log_info("Saved vtk file as %s\n", save_path);
                }

                return;
            } else if(kp == KEY_F) {
                gui_state->search_window.show = true;
                gui_state->search_window.bounds.x = (float)GetScreenWidth() / 2.0f - gui_state->search_window.bounds.width;
                gui_state->search_window.bounds.y = (float)GetScreenHeight() / 2.0f - gui_state->search_window.bounds.height;
                return;
            }
        }

        if(kp >= KEY_ZERO && kp <= KEY_NINE) {
            int index = kp - KEY_ONE - 1;

            if(index == -2) {
                // zero pressed. we consider 8 (10th column)
                index = 8;
            }

            if(index <= gui_state->max_data_index - 1) {
                gui_state->current_data_index = index;
            }
        }

        if(kp == KEY_PAGE_UP) {
            int index = gui_state->current_data_index + 1;

            if(index <= gui_state->max_data_index - 1) {
                gui_state->current_data_index = index;
            }

        } else if(kp == KEY_PAGE_DOWN) {
            int index = gui_state->current_data_index - 1;

            if(index >= -1) {
                gui_state->current_data_index = index;
            }
        }

        if(kp == KEY_S) {
            if(IN_DRAW && gui_config->enable_slice && gui_state->current_mode != EDITING) {
                gui_state->current_mode = SLICING;
                gui_state->slicing_mesh = false;

                if(gui_state->old_cell_visibility) {
                    reset_grid_visibility(gui_config, gui_state);
                }
            }

            if(IN_DRAW && gui_state->current_mode == EDITING) {
                const char *save_path;

                if(gui_state->ctrl_pressed) {
                    char const *filter[1] = {"*.alg"};
                    save_path = tinyfd_saveFileDialog("Save ALG file", gui_config->input, 1, filter, "alg files");
                } else {
                    save_path = gui_config->grid_info.file_name;
                }

                if(save_path) {
                    save_vtk_unstructured_grid_as_alg_file(gui_config->grid_info.vtk_grid, save_path, false);
                    log_info("Saved vtk file as %s\n", save_path);
                }
            }
        }

        if(gui_state->current_mode == SLICING) {

            if(IsKeyDown(KEY_ENTER)) {
                gui_state->current_mode = VISUALIZING;
                gui_state->slicing_mesh = true;
            }

            if(IsKeyDown(KEY_BACKSPACE)) {
                gui_state->current_mode = VISUALIZING;
                gui_state->slicing_mesh = false;
                reset_grid_visibility(gui_config, gui_state);
            }

            float slice_inc = 0.8f;

            if(IsKeyDown(KEY_LEFT_ALT)) {
                // Plane roll (z-axis) controls
                if(IsKeyDown(KEY_LEFT))
                    gui_state->plane_roll += slice_inc;
                else if(IsKeyDown(KEY_RIGHT))
                    gui_state->plane_roll -= slice_inc;

                if(IsKeyDown(KEY_DOWN))
                    gui_state->plane_pitch += slice_inc;
                else if(IsKeyDown(KEY_UP))
                    gui_state->plane_pitch -= slice_inc;

            } else {

                slice_inc = 0.05f;

                if(IsKeyDown(KEY_LEFT_CONTROL)) {
                    slice_inc = 0.01f;
                }
                // Plane roll (z-axis) controls
                if(IsKeyDown(KEY_LEFT))
                    gui_state->plane_tx -= slice_inc;
                else if(IsKeyDown(KEY_RIGHT))
                    gui_state->plane_tx += slice_inc;

                if(IsKeyDown(KEY_DOWN))
                    gui_state->plane_ty -= slice_inc;
                else if(IsKeyDown(KEY_UP))
                    gui_state->plane_ty += slice_inc;
            }
        } else {
            if(kp == KEY_RIGHT || IsKeyDown(KEY_UP)) {

                if(gui_config->draw_type == DRAW_FILE && gui_config->final_file_index == 0)
                    return;

                gui_config->current_file_index++;
                CHECK_FILE_INDEX(gui_config)
                omp_unset_nest_lock(&gui_config->sleep_lock);
                return;
            }

            if(gui_config->draw_type == DRAW_FILE) {
                // Return one step only works on file visualization...
                if(kp == KEY_LEFT || IsKeyDown(KEY_DOWN)) {
                    if(gui_config->final_file_index == 0)
                        return;
                    gui_config->current_file_index--;
                    CHECK_FILE_INDEX(gui_config)
                    omp_unset_nest_lock(&gui_config->sleep_lock);
                    return;
                }
            }
        }
    }

    if(kp == KEY_Q) {
        gui_state->scale.window.show = !gui_state->scale.window.show;
        return;
    }

    if(IsKeyDown(KEY_A)) {
        if(gui_state->voxel_alpha - 1 >= 1) {
            gui_state->voxel_alpha = gui_state->voxel_alpha - 1;
        }
        return;
    }

    if(IsKeyDown(KEY_Z)) {
        if(gui_state->voxel_alpha + 1 <= 255) {
            gui_state->voxel_alpha = gui_state->voxel_alpha + 1;
        }
        return;
    }

    if(kp == KEY_G) {
        gui_state->draw_grid_only = !gui_state->draw_grid_only;
        return;
    }

    if(kp == KEY_PERIOD) {
        gui_state->current_scale = (gui_state->current_scale + 1) % NUM_SCALES;
        return;
    }

    if(kp == KEY_COMMA) {
        if(gui_state->current_scale - 1 >= 0) {
            gui_state->current_scale = (gui_state->current_scale - 1);
        } else {
            gui_state->current_scale = NUM_SCALES - 1;
        }
        return;
    }

    if(kp == KEY_L) {
        gui_state->draw_grid_lines = !gui_state->draw_grid_lines;
        return;
    }

    if(kp == KEY_SPACE) {
        if(gui_state->current_mode == VISUALIZING) {
            if(gui_config->draw_type == DRAW_FILE) {
                if(gui_config->final_file_index != 0) {
                    gui_config->paused = !gui_config->paused;
                }
            } else {
                gui_config->paused = !gui_config->paused;
            }
            return;
        }
    }

    if(kp == KEY_R) {
        bool full_reset = false;

        if(IsKeyDown(KEY_LEFT_ALT)) {
            full_reset = true;
        }

        reset(gui_config, gui_state, full_reset);
        return;
    }

    if(kp == KEY_X) {
        gui_state->ap_graph_config->graph.show = !gui_state->ap_graph_config->graph.show;
        return;
    }

    if(kp == KEY_C) {
        gui_state->scale.window.show = gui_state->c_pressed;
        gui_state->ap_graph_config->graph.show = gui_state->c_pressed;
        gui_state->end_info_box.window.show = gui_state->c_pressed;
        gui_state->mesh_info_box.window.show = gui_state->c_pressed;
        gui_state->controls_window.show = gui_state->c_pressed;

        gui_state->c_pressed = !gui_state->c_pressed;
        gui_state->show_coordinates = !gui_state->show_coordinates;
        return;
    }

    if(kp == KEY_H) {
        if(gui_state->current_mode == VISUALIZING) {
            gui_state->help_box.window.show = !gui_state->help_box.window.show;
        }
        if(gui_state->current_mode == SLICING) {
            gui_state->slice_help_box.window.show = !gui_state->slice_help_box.window.show;
        } else if(gui_state->current_mode == EDITING) {
            gui_state->edit_help_box.window.show = !gui_state->edit_help_box.window.show;
        }
        return;
    }

    if(gui_config->draw_type == DRAW_FILE) {
        if(kp == KEY_O) {
            if(!gui_config->paused)
                return;

            char *buf = get_current_directory();

            char const *tmp = tinyfd_selectFolderDialog("Select a directory", buf);
            if(tmp) {
                gui_config->input = strdup(tmp);
                reset(gui_config, gui_state, true);
                free(gui_config->message);
                gui_config->message = strdup("Loading Mesh...");

            } else {
                gui_config->input = NULL;
            }

            free(buf);

            return;
        }

        if(kp == KEY_F) {

            gui_config->paused = true;

            char *buf = get_current_directory();

            char const *filter[] = {"*.geo", "*.Esca", "*.pvd", "*.acm", "*.vtk", "*.vtu"};

            char const *tmp = tinyfd_openFileDialog("Select a simulation file", buf, 4, filter, "simulation result (pvd, vtk, vtu or acm)", 0);

            if(tmp) {
                gui_config->input = strdup(tmp);
            } else {
                gui_config->input = NULL;
            }

            free(buf);

            if(tmp) {
                reset(gui_config, gui_state, true);
                free(gui_config->message);
                gui_config->message = strdup("Loading Mesh...");
            }
            return;
        }
    }
}

static void handle_input(struct gui_shared_info *gui_config, struct gui_state *gui_state) {

    gui_state->mouse_pos = GetMousePosition();

    if(gui_state->handle_keyboard_input) {
        handle_keyboard_input(gui_config, gui_state);
    }

    if(IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {

        gui_state->ray = GetMouseRay(gui_state->mouse_pos, gui_state->camera);

        if(gui_state->current_mode == EDITING) {
            gui_state->get_cell_property = true;
        } else {
            if(!gui_state->search_window.show) {
                if(gui_state->mouse_timer == -1) {
                    gui_state->double_clicked = false;
                    gui_state->mouse_timer = GetTime();
                } else {
                    double delay = GetTime() - gui_state->mouse_timer;
                    if(delay < DOUBLE_CLICK_DELAY) {
                        gui_state->double_clicked = true;
                        gui_state->mouse_timer = -1;
                    } else {
                        gui_state->mouse_timer = -1;
                        gui_state->double_clicked = false;
                    }
                }
            }
        }

        check_colisions_for_move(gui_state);

    } else if(IsMouseButtonPressed(MOUSE_RIGHT_BUTTON)) {

        if(gui_state->current_mode == EDITING) {
            gui_state->ray = GetMouseRay(gui_state->mouse_pos, gui_state->camera);
            gui_state->paste_cell_property = true;
        } else {
            if(hmlen(gui_state->ap_graph_config->selected_aps) && gui_state->ap_graph_config->graph.show) {
                if(CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->graph.bounds)) {
                    if(gui_state->ap_graph_config->selected_point_for_apd1.x == FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd1.y == FLT_MAX) {
                        gui_state->ap_graph_config->selected_point_for_apd1.x = gui_state->mouse_pos.x;
                        gui_state->ap_graph_config->selected_point_for_apd1.y = gui_state->mouse_pos.y;
                    } else {
                        if(gui_state->ap_graph_config->selected_point_for_apd2.x == FLT_MAX &&
                           gui_state->ap_graph_config->selected_point_for_apd2.y == FLT_MAX) {
                            gui_state->ap_graph_config->selected_point_for_apd2.x = gui_state->mouse_pos.x;
                            gui_state->ap_graph_config->selected_point_for_apd2.y = gui_state->ap_graph_config->selected_point_for_apd1.y;
                        } else {
                            gui_state->ap_graph_config->selected_point_for_apd1.x = gui_state->mouse_pos.x;
                            gui_state->ap_graph_config->selected_point_for_apd1.y = gui_state->mouse_pos.y;

                            gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
                            gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;
                        }
                    }
                }
            }
        }
    }

    maybe_move_or_drag(gui_state);

    if(hmlen(gui_state->ap_graph_config->selected_aps) && gui_state->ap_graph_config->graph.show) {
        float t = Remap(gui_state->mouse_pos.x, gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, 0.0f, gui_config->final_time);
        float v = Remap(gui_state->mouse_pos.y, gui_state->ap_graph_config->min_y, gui_state->ap_graph_config->max_y, gui_config->min_v, gui_config->max_v);
        if(CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->graph.bounds)) {
            gui_state->ap_graph_config->selected_ap_point.x = t;
            gui_state->ap_graph_config->selected_ap_point.y = v;
        } else {
            gui_state->ap_graph_config->selected_ap_point.x = FLT_MAX;
            gui_state->ap_graph_config->selected_ap_point.y = FLT_MAX;
        }
    }
}

static void configure_info_boxes_sizes(struct gui_state *gui_state, int help_box_lines, int slice_help_box_lines, int edit_help_box_lines,
                                       int end_info_box_lines, float box_w, float text_offset) {

    gui_state->help_box.window.bounds.width = box_w;
    gui_state->help_box.window.bounds.height = (text_offset * (float)help_box_lines) + BOX_MARGIN + 10.0f;
    gui_state->help_box.num_lines = help_box_lines;

    gui_state->slice_help_box.window.bounds.width = box_w - 100;
    gui_state->slice_help_box.window.bounds.height = (text_offset * (float)slice_help_box_lines) + BOX_MARGIN + 10.0f;
    gui_state->slice_help_box.num_lines = slice_help_box_lines;

    gui_state->edit_help_box.window.bounds.width = box_w - 60;
    gui_state->edit_help_box.window.bounds.height = (text_offset * (float)edit_help_box_lines) + BOX_MARGIN + 10.0f;
    gui_state->edit_help_box.num_lines = edit_help_box_lines;

    box_w = box_w - 100;
    gui_state->mesh_info_box.window.bounds.width = box_w;

    gui_state->mesh_info_box.window.bounds.x = (float)gui_state->current_window_width - box_w - BOX_MARGIN;
    gui_state->mesh_info_box.window.bounds.y = 10.0f;

    gui_state->end_info_box.window.bounds.width = box_w;
    gui_state->end_info_box.window.bounds.height = (text_offset * (float)end_info_box_lines) + BOX_MARGIN + 10.0f;
    gui_state->end_info_box.num_lines = end_info_box_lines;

    gui_state->end_info_box.window.bounds.x = gui_state->mesh_info_box.window.bounds.x - gui_state->mesh_info_box.window.bounds.width - BOX_MARGIN;
    gui_state->end_info_box.window.bounds.y = gui_state->mesh_info_box.window.bounds.y;
}

void init_and_open_gui_window(struct gui_shared_info *gui_config) {

    const int end_info_box_lines = 9;
    const int font_size_small = 16;
    const int font_size_big = 20;
    const int help_box_lines = SIZEOF(help_box_strings);
    const int slice_help_box_lines = SIZEOF(slice_help_box_strings);
    const int edit_help_box_lines = SIZEOF(edit_help_box_strings);

    omp_set_nest_lock(&gui_config->sleep_lock);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT);
    char window_title[4096];

    const int draw_type = gui_config->draw_type;

    if(draw_type == DRAW_SIMULATION) {
        sprintf(window_title, "Simulation visualization - %s", gui_config->config_name);
    } else {
        sprintf(window_title, "Monoalg3D Simulator");
    }

    InitWindow(0, 0, window_title);

    SetTargetFPS(120);

    Image icon = LoadImage("res/icon.png");

    if(icon.data) {
        SetWindowIcon(icon);
    }

    UnloadImage(icon);

    if(gui_config->ui_scale == 0.0) {
        Vector2 ui_scale = GetWindowScaleDPI();
        gui_config->ui_scale = ui_scale.x;
    }

    MaximizeWindow();

    struct gui_state *gui_state = new_gui_state_with_font_sizes((float)font_size_small, (float)font_size_big, gui_config->ui_scale);

    gui_state->help_box.lines = (char **)help_box_strings;
    gui_state->help_box.title = "Default controls";

    gui_state->slice_help_box.lines = (char **)slice_help_box_strings;
    gui_state->slice_help_box.title = "Mesh slicing help";

    gui_state->edit_help_box.lines = (char **)edit_help_box_strings;
    gui_state->edit_help_box.title = "Mesh slicing help";

    float wider_text_w = 0;
    float wider_text_h = 0;

    for(int i = 0; i < help_box_lines; i++) {
        Vector2 tmp = MeasureTextEx(gui_state->font, help_box_strings[i], gui_state->font_size_small, gui_state->font_spacing_small);
        if(tmp.x > wider_text_w) {
            wider_text_w = tmp.x;
            wider_text_h = tmp.y;
        }
    }

    const float text_offset = 1.5f * wider_text_h;
    const float box_w = wider_text_w * 1.08f;

    gui_state->end_info_box.lines = (char **)malloc(sizeof(char *) * end_info_box_lines);
    gui_state->end_info_box.title = "Simulation time";

    gui_state->mesh_info_box.title = "Mesh information";
    gui_state->mesh_info_box.lines = NULL;
    configure_info_boxes_sizes(gui_state, help_box_lines, slice_help_box_lines, edit_help_box_lines, end_info_box_lines, box_w, text_offset);

    Vector2 message_width;
    struct mesh_info *mesh_info = new_mesh_info();
    bool end_info_box_strings_configured = false;

    Model plane;
    Color plane_color = WHITE;
    plane_color.a = 100;

    if(NOT_IN_DRAW) {
        gui_config->progress = 100;
    }

    void (*draw_grid_function)(struct gui_shared_info *, struct gui_state *, int, struct draw_context *);
    void (*draw_purkinje_function)(struct gui_shared_info *, struct gui_state *);

    draw_purkinje_function = NULL;

    if(draw_type == DRAW_SIMULATION) {
        draw_grid_function = &draw_alg_mesh;
        draw_purkinje_function = &draw_alg_purkinje_network;
    } else {
        draw_grid_function = &draw_vtk_unstructured_grid;
        draw_purkinje_function = &draw_vtk_purkinje_network;
    }

    struct draw_context draw_context = {0};

    draw_context.shader = LoadShader("res/instanced_vertex_shader.vs", "res/fragment_shader.fs");
    draw_context.shader.locs[SHADER_LOC_MATRIX_MVP] = GetShaderLocation(draw_context.shader, "mvp");
    draw_context.shader.locs[SHADER_LOC_MATRIX_MODEL] = GetShaderLocationAttrib(draw_context.shader, "instanceTransform");
    draw_context.shader.locs[SHADER_LOC_VERTEX_COLOR] = GetShaderLocationAttrib(draw_context.shader, "color");
    draw_context.shader.locs[SHADER_LOC_VECTOR_VIEW] = GetShaderLocation(draw_context.shader, "viewPos");
    draw_context.mesh = GenMeshCube(1.0f, 1.0f, 1.0f);

    // Lights
    const float light_offset = 0.6f;
    Vector3 light_pos = gui_state->camera.position;
    light_pos.z = light_pos.z + light_offset;
    gui_state->light = CreateLight(LIGHT_DIRECTIONAL, light_pos, gui_state->camera.target, WHITE, draw_context.shader);

    Matrix rotation_matrix;

    int grid_mask;

    while(!WindowShouldClose()) {

        if(gui_state->draw_grid_only) {
            grid_mask = 2;
        } else if(gui_state->draw_grid_lines) {
            grid_mask = 1;
        } else {
            grid_mask = 0;
        }

        UpdateCamera(&(gui_state->camera));

        light_pos = gui_state->camera.position;
        light_pos.z = light_pos.z + light_offset;
        gui_state->light.position = light_pos;
        gui_state->light.target = gui_state->camera.target;
        UpdateLightValues(draw_context.shader, gui_state->light);

        BeginDrawing();

        if(IsWindowResized()) {
            gui_state->current_window_width = GetScreenWidth();
            gui_state->current_window_height = GetScreenHeight();
            reset_ui(gui_state);
        }

        handle_input(gui_config, gui_state);

        if(gui_config->grid_info.loaded) {

            omp_set_nest_lock(&gui_config->draw_lock);

            uint32_t n_active;

            if(draw_type == DRAW_SIMULATION) {
                n_active = gui_config->grid_info.alg_grid->num_active_cells;
            } else {
                n_active = gui_config->grid_info.vtk_grid->num_cells;
            }

            bool realloc_matrices = gui_state->last_n_active != n_active;
            gui_state->last_n_active = n_active;

            if(draw_type == DRAW_FILE) {
                if(gui_config->grid_info.file_name) {
                    sprintf(window_title, "Visualizing file - %s", gui_config->grid_info.file_name);
                    SetWindowTitle(window_title);
                }
            }
            ClearBackground(GRAY);

            if(!mesh_info->center_calculated) {
                if(draw_type == DRAW_SIMULATION) {
                    gui_state->mesh_offset = find_mesh_center(gui_config->grid_info.alg_grid, mesh_info);
                    gui_state->max_data_index = 0;
                } else {
                    gui_state->mesh_offset = find_mesh_center_vtk(gui_config->grid_info.vtk_grid, mesh_info);
                    gui_state->max_data_index = arrlen(gui_config->grid_info.vtk_grid->extra_values);
                }

                gui_state->mesh_scale_factor = fmaxf(gui_state->mesh_offset.x, fmaxf(gui_state->mesh_offset.y, gui_state->mesh_offset.z)) / 1.8f;

                const float scale = gui_state->mesh_scale_factor;

                Vector2 pos;
                pos.x = (mesh_info->max_size.x - mesh_info->min_size.x) / scale;
                pos.y = (mesh_info->max_size.z - mesh_info->min_size.z) / scale;

                const float mult = 1.2f;
                plane = LoadModelFromMesh(GenMeshCube(pos.x * mult, 0.1f / scale, pos.y * mult));
            }

            BeginMode3D(gui_state->camera);

            if(gui_state->current_mode == SLICING) {
                plane.transform = MatrixTranslate(gui_state->plane_tx, gui_state->plane_ty, gui_state->plane_tz);
                rotation_matrix = MatrixRotateXYZ((Vector3){DEG2RAD * gui_state->plane_pitch, 0.0f, DEG2RAD * gui_state->plane_roll});
                plane.transform = MatrixMultiply(plane.transform, rotation_matrix);

                gui_state->plane_normal = Vector3Normalize(Vector3Transform((Vector3){0, 1, 0}, rotation_matrix));
                gui_state->plane_point = Vector3Transform((Vector3){0, 0, 0}, plane.transform);

                DrawModel(plane, (Vector3){0, 0, 0}, 1.0f, plane_color);
            }

            if(!gui_state->slicing_mesh) {

                if(gui_config->grid_info.loaded && realloc_matrices) {

                    free(draw_context.translations);
                    free(draw_context.colors);
                    free(draw_context.instance_transforms);
                    free(draw_context.colors_transforms);

                    draw_context.colors = malloc(n_active * sizeof(Color));
                    draw_context.translations = malloc(n_active * sizeof(Matrix)); // Locations of instances

                    draw_context.instance_transforms = (float16 *)malloc(n_active * sizeof(float16));
                    draw_context.colors_transforms = (float4 *)malloc(n_active * sizeof(float4));

                    gui_state->handle_keyboard_input = true;
                }

                draw_grid_function(gui_config, gui_state, grid_mask, &draw_context);

                if(draw_purkinje_function != NULL) {
                    draw_purkinje_function(gui_config, gui_state);
                }
            }

            gui_state->double_clicked = false;

            if(gui_state->show_coordinates) {
                draw_coordinates(gui_state);
            }

            EndMode3D();

            if(gui_state->slicing_mesh) {
                float spacing = gui_state->font_spacing_big;
                static Color c = RED;

                message_width = MeasureTextEx(gui_state->font, "Slicing Mesh...", gui_state->font_size_big, spacing);
                int posx = GetScreenWidth() / 2 - (int)message_width.x / 2;
                int posy = GetScreenHeight() / 2 - 50;

                int rec_width = (int)(message_width.x) + 50;
                int rec_height = (int)(message_width.y) + 2;

                DrawRectangle(posx, posy, rec_width, rec_height, c);

                DrawTextEx(gui_state->font, "Slicing Mesh...", (Vector2){(float)posx + ((float)rec_width - message_width.x) / 2, (float)posy},
                           gui_state->font_size_big, gui_state->font_spacing_big, BLACK);
            }

            // We finished drawing everything that depends on the mesh being loaded
            omp_unset_nest_lock(&gui_config->draw_lock);

            if(gui_state->controls_window.show) {
                draw_control_window(gui_state, gui_config);
            }

            if(gui_state->show_coordinates) {
                DrawText("x", (int)gui_state->coordinates_label_x_position.x, (int)gui_state->coordinates_label_x_position.y, (int)gui_state->font_size_big,
                         RED);
                DrawText("y", (int)gui_state->coordinates_label_y_position.x, (int)gui_state->coordinates_label_y_position.y, (int)gui_state->font_size_big,
                         GREEN);
                DrawText("z", (int)gui_state->coordinates_label_z_position.x, (int)gui_state->coordinates_label_z_position.y, (int)gui_state->font_size_big,
                         DARKBLUE);
            }

            if(gui_state->mesh_info_box.window.show) {
                bool configured = configure_mesh_info_box_strings(gui_state, gui_config, &gui_state->mesh_info_box.lines, draw_type, mesh_info);

                if(configured) {
                    int mesh_info_box_lines = arrlen(gui_state->mesh_info_box.lines);
                    gui_state->mesh_info_box.num_lines = mesh_info_box_lines;
                    gui_state->mesh_info_box.window.bounds.height = (text_offset * (float)mesh_info_box_lines) + BOX_MARGIN + 10.0f;
                    draw_text_window(&gui_state->mesh_info_box, gui_state, text_offset);

                    for(int i = 0; i < mesh_info_box_lines; i++) {
                        free(gui_state->mesh_info_box.lines[i]);
                    }
                    arrfree(gui_state->mesh_info_box.lines);
                    gui_state->mesh_info_box.lines = NULL;
                }
            }

            if(gui_state->scale.window.show) {
                draw_scale(gui_config->min_v, gui_config->max_v, gui_state, gui_config->int_scale);
            }

            if(gui_state->ap_graph_config->graph.show && hmlen(gui_state->ap_graph_config->selected_aps)) {
                draw_ap_graph(gui_state, gui_config);
            }

            if(gui_state->help_box.window.show) {
                draw_text_window(&gui_state->help_box, gui_state, text_offset);
            }

            if(gui_state->slice_help_box.window.show) {
                draw_text_window(&gui_state->slice_help_box, gui_state, text_offset);
            }

            if(gui_state->edit_help_box.window.show) {
                draw_text_window(&gui_state->edit_help_box, gui_state, text_offset);
            }

            if(!gui_config->simulating) {
                if(draw_type == DRAW_SIMULATION) {
                    if(gui_state->end_info_box.window.show) {
                        if(!end_info_box_strings_configured) {
                            configure_end_info_box_strings(gui_config, &gui_state->end_info_box.lines);
                            end_info_box_strings_configured = true;
                        }
                        draw_text_window(&gui_state->end_info_box, gui_state, text_offset);
                    }
                }
            }

            if(!gui_config->paused) {
                omp_unset_nest_lock(&gui_config->sleep_lock);
            }

            if(gui_state->search_window.show) {
                bool update = !draw_search_window(gui_state);
                gui_state->search_window.show = update;
                gui_state->handle_keyboard_input = !update;
            }
        } else {

            gui_state->handle_keyboard_input = false;
            ClearBackground(GRAY);
            float spacing = gui_state->font_spacing_big;
            static Color c = RED;

            if(!gui_config->message) {
                gui_config->message = strdup("Loading Mesh...");
                c = WHITE;
            }

            message_width = MeasureTextEx(gui_state->font, gui_config->message, gui_state->font_size_big, spacing);
            int posx = GetScreenWidth() / 2 - (int)message_width.x / 2;
            int posy = GetScreenHeight() / 2 - 50;

            int rec_width = (int)(message_width.x) + 50;
            int rec_height = (int)(message_width.y) + 2;

            int rec_bar_w = (int)Remap((float)gui_config->progress, 0, (float)gui_config->file_size, 0, (float)rec_width);

            DrawRectangle(posx, posy, rec_bar_w, rec_height, c);
            DrawRectangleLines(posx, posy, rec_width, rec_height, BLACK);

            // This should not happen... but it does....
            if(gui_config->message) {
                DrawTextEx(gui_state->font, gui_config->message, (Vector2){(float)posx + ((float)rec_width - message_width.x) / 2, (float)posy},
                           gui_state->font_size_big, gui_state->font_spacing_big, BLACK);

                gui_state->handle_keyboard_input = true;
            }
        }

        // Draw FPS
        int fps = GetFPS();
        const char *text = TextFormat("%2i FPS - Frame Time %lf", fps, GetFrameTime());
        Vector2 text_size = MeasureTextEx(gui_state->font, text, gui_state->font_size_big, gui_state->font_spacing_big);

        DrawTextEx(gui_state->font, text,
                   (Vector2){((float)gui_state->current_window_width - text_size.x - 10.0f), ((float)gui_state->current_window_height - text_size.y - 30)},
                   gui_state->font_size_big, gui_state->font_spacing_big, BLACK);

        text_size = MeasureTextEx(gui_state->font, "Press H to show/hide the help box", gui_state->font_size_big, gui_state->font_spacing_big);

        if(gui_state->recalculating_visibility) {
            DrawTextEx(gui_state->font, "Recalculating mesh visibility", (Vector2){10.0f, ((float)gui_state->current_window_height - text_size.y - 30.0f)},
                       gui_state->font_size_big, gui_state->font_spacing_big, BLACK);

        } else {
            if(gui_state->current_mode == VISUALIZING) {
                DrawTextEx(gui_state->font, "Press H to show/hide the help box",
                           (Vector2){10.0f, ((float)gui_state->current_window_height - text_size.y - 30.0f)}, gui_state->font_size_big,
                           gui_state->font_spacing_big, BLACK);
            } else if(gui_state->current_mode == SLICING) {
                DrawTextEx(gui_state->font, "Slicing mode - Press H to show/hide the help box.\nPress Enter to confirm or Backspace to reset and exit.",
                           (Vector2){10.0f, ((float)gui_state->current_window_height - text_size.y - 35.0f)}, gui_state->font_size_big,
                           gui_state->font_spacing_big, BLACK);
            } else if(gui_state->current_mode == EDITING) {
                DrawTextEx(gui_state->font,
                           "Editing mode - Press H to show/hide the help box.\nPress E exit edit mode or S to save the changes to the mesh file.",
                           (Vector2){10.0f, ((float)gui_state->current_window_height - text_size.y - 35.0f)}, gui_state->font_size_big,
                           gui_state->font_spacing_big, BLACK);
            }
        }

        float upper_y = text_size.y + 30;

        if(gui_state->ctrl_pressed && gui_state->current_mouse_over_volume.position_draw.x != -1) {

            Vector2 info_pos;

            if(gui_config->draw_type == DRAW_SIMULATION) {
                text = TextFormat("Mouse is on Volume: %.2lf, %.2lf, %.2lf with grid position %i and value %.2lf",
                                  gui_state->current_mouse_over_volume.position_draw.x, gui_state->current_mouse_over_volume.position_draw.y,
                                  gui_state->current_mouse_over_volume.position_draw.z,
                                  gui_state->current_mouse_over_volume.matrix_position + gui_state->current_mouse_over_volume.v);

            } else {

                text = TextFormat("Mouse is on Volume: %.2lf, %.2lf, %.2lf with value %.2lf", gui_state->current_mouse_over_volume.position_draw.x,
                                  gui_state->current_mouse_over_volume.position_draw.y, gui_state->current_mouse_over_volume.position_draw.z,
                                  gui_state->current_mouse_over_volume.v);
            }

            text_size = MeasureTextEx(gui_state->font, text, gui_state->font_size_big, gui_state->font_spacing_big);

            info_pos =
                (Vector2){((float)gui_state->current_window_width - text_size.x - 10), ((float)gui_state->current_window_height - text_size.y - upper_y)};

            DrawTextEx(gui_state->font, text, info_pos, gui_state->font_size_big, gui_state->font_spacing_big, BLACK);
        }

        EndDrawing();

        // TODO: this should go to a separate thread (maybe??)
        if(gui_state->slicing_mesh) {
            set_visibility_after_split(gui_config, gui_state);
            gui_state->slicing_mesh = false;
        }
    }

    gui_config->exit = true;
    free(mesh_info);

    if(end_info_box_strings_configured) {
        for(int i = 0; i < end_info_box_lines; i++) {
            free(gui_state->end_info_box.lines[i]);
        }
    }

    free(gui_state->end_info_box.lines);

    hmfree(gui_state->ap_graph_config->selected_aps);

    free(gui_state->ap_graph_config);
    free(gui_state);

    omp_unset_nest_lock(&gui_config->draw_lock);
    omp_unset_nest_lock(&gui_config->sleep_lock);

    free(draw_context.translations);
    free(draw_context.colors);
    free(draw_context.instance_transforms);
    free(draw_context.colors_transforms);

    CloseWindow();
}
