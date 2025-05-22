#include "gui_draw.h"
#include "gui.h"
#include <float.h>

void set_camera_params(Camera3D *camera, bool set_mode) {
    camera->position = (Vector3){0.1f, 0.1f, 12.0f};
    camera->target = (Vector3){0.f, 0.f, 0.f};
    camera->up = (Vector3){0.0f, 1.0f, 0.0f};
    camera->fovy = 45.0f;

    if(set_mode) {
        camera->projection = CAMERA_PERSPECTIVE;
        SetCameraMode(*camera, CAMERA_FREE);
    }
}

void draw_ap_graph(struct gui_state *gui_state, struct gui_shared_info *gui_config) {

    int n = hmlen(gui_state->ap_graph_config->selected_aps);

    const float graph_w = gui_state->ap_graph_config->graph.bounds.width;
    const float graph_h = gui_state->ap_graph_config->graph.bounds.height;

    if(gui_state->ap_graph_config->graph.bounds.x + graph_w > (float)gui_state->current_window_width || gui_state->ap_graph_config->graph.bounds.x < 0) {
        gui_state->ap_graph_config->graph.bounds.x = 10;
    }

    if(gui_state->ap_graph_config->graph.bounds.y + graph_h > (float)gui_state->current_window_height) {
        gui_state->ap_graph_config->graph.bounds.y = (float)gui_state->current_window_height - graph_h - 90.0f;
    }

    if(gui_state->ap_graph_config->graph.bounds.y < 0) {
        gui_state->ap_graph_config->graph.bounds.y = 0;
    }

    const float graph_x = gui_state->ap_graph_config->graph.bounds.x;
    const float graph_y = gui_state->ap_graph_config->graph.bounds.y;

    static const Color colors[] = {DARKGRAY, GOLD,     ORANGE, PINK,   RED,        MAROON, GREEN,     LIME,  DARKGREEN,
                                   BLUE,     DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BROWN,  DARKBROWN, BLACK, MAGENTA};

    int num_colors = SIZEOF(colors);

    Font font = gui_state->font;

    float font_size_big = gui_state->font_size_big - 6;
    float font_size_small = gui_state->font_size_small - 6;

    float spacing_big = font_size_big / (float)font.baseSize;
    float spacing_small = font_size_small / (float)font.baseSize;

    Rectangle graph_window = gui_state->ap_graph_config->graph.bounds;
    graph_window.y -= WINDOW_STATUSBAR_HEIGHT;
    graph_window.height += WINDOW_STATUSBAR_HEIGHT;
    gui_state->ap_graph_config->graph.show = !GuiWindowBox(graph_window, TextFormat("Selected APs (%i)", n));

    gui_state->ap_graph_config->drag_graph_button = (Rectangle){(graph_x + graph_w) - 16.0f, graph_y + graph_h - 16.0f, 16.0f, 16.0f};

    GuiButton(gui_state->ap_graph_config->drag_graph_button, "#71#");

    char tmp[TMP_SIZE];

    snprintf(tmp, TMP_SIZE, "%.2lf", gui_config->final_time);
    Vector2 text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    snprintf(tmp, TMP_SIZE, "%.2lf", gui_config->min_v);
    Vector2 text_width_y = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    gui_state->ap_graph_config->min_x = graph_x + 2.6f * text_width_y.x;
    gui_state->ap_graph_config->max_x = graph_x + graph_w - text_width.x;

    gui_state->ap_graph_config->min_y = graph_y + graph_h - text_width.y;
    gui_state->ap_graph_config->max_y = graph_y + 20.0f; // This is actually the smallest allowed y

    Vector2 p1, p2;

    uint32_t num_ticks;
    float tick_ofsset = 10;
    num_ticks = (int)(gui_config->final_time / tick_ofsset);

    if(num_ticks < MIN_HORIZONTAL_TICKS) {
        num_ticks = MIN_HORIZONTAL_TICKS;
        tick_ofsset = gui_config->final_time / (float)num_ticks;
    } else if(num_ticks > MAX_HORIZONTAL_TICKS) {
        num_ticks = MAX_HORIZONTAL_TICKS;
        tick_ofsset = gui_config->final_time / (float)num_ticks;
    }

    char *time_text;
    char *time_template;
    bool steps = false;

    if(gui_config->dt == 0) {
        time_text = "Time (steps)";
        time_template = "%d";
        steps = true;
    } else {
        time_text = "Time (ms)";
        time_template = "%.2lf";
    }

    // Draw x label
    text_width = MeasureTextEx(font, time_text, font_size_big, spacing_big);

    gui_state->ap_graph_config->min_y -= text_width.y * 1.5f;

    float min_x = gui_state->ap_graph_config->min_x;
    float min_y = gui_state->ap_graph_config->min_y;

    float max_x = gui_state->ap_graph_config->max_x;
    float max_y = gui_state->ap_graph_config->max_y;

    Vector2 text_position = (Vector2){graph_x + gui_state->ap_graph_config->graph.bounds.width / 2.0f - text_width.x / 2.0f, min_y + text_width.y};

    DrawTextEx(font, time_text, text_position, font_size_big, spacing_big, BLACK);

    float time = 0.0f;

    // Draw horizontal ticks (t)
    for(uint32_t t = 0; t <= num_ticks; t++) {

        p1.x = Remap(time, 0.0f, gui_config->final_time, min_x, max_x);
        p1.y = min_y - 5;

        p2.x = p1.x;
        p2.y = min_y + 5;

        if(!(t % 2)) {

            float remaped_data = Remap(p1.x, min_x, max_x, 0.0f, gui_config->final_time);

            if(steps) {
                snprintf(tmp, TMP_SIZE, time_template, (int)remaped_data);
            } else {
                snprintf(tmp, TMP_SIZE, time_template, remaped_data);
            }

            text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);
            DrawTextEx(font, tmp, (Vector2){p1.x - text_width.x / 2.0f, p1.y + 10}, font_size_small, spacing_small, RED);
        }

        DrawLineV(p1, p2, RED);

        DrawLineV(p1, (Vector2){p1.x, gui_state->ap_graph_config->max_y}, LIGHTGRAY);
        time += tick_ofsset;
    }

    tick_ofsset = 10;
    num_ticks = (uint32_t)((gui_config->max_v - gui_config->min_v) / tick_ofsset);

    if(num_ticks < MIN_VERTICAL_TICKS) {
        num_ticks = MIN_VERTICAL_TICKS;
        tick_ofsset = (gui_config->max_v - gui_config->min_v) / (float)num_ticks;
    } else if(num_ticks > MAX_VERTICAL_TICKS) {
        num_ticks = MAX_VERTICAL_TICKS;
        tick_ofsset = (gui_config->max_v - gui_config->min_v) / (float)num_ticks;
    }

    // Draw y label
    {
        char *label;
        label = "Vm (mV)";
        text_width = MeasureTextEx(font, label, font_size_big, spacing_big);

        text_position.x = gui_state->ap_graph_config->graph.bounds.x + 15;
        text_position.y = min_y - ((min_y - max_y) / 2.0f) + (text_width.x / 2.0f);

        DrawTextPro(font, label, text_position, (Vector2){0, 0}, -90, font_size_big, spacing_big, BLACK);
    }

    float v = (float)gui_config->min_v;
    snprintf(tmp, TMP_SIZE, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    // Draw vertical ticks (Vm)
    for(uint32_t t = 0; t <= num_ticks; t++) {

        p1.x = (float)graph_x + 5.0f;
        p1.y = Remap(v, gui_config->min_v, gui_config->max_v, min_y, gui_state->ap_graph_config->max_y);

        snprintf(tmp, TMP_SIZE, "%.2lf", Remap(p1.y, min_y, gui_state->ap_graph_config->max_y, gui_config->min_v, gui_config->max_v));
        text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

        DrawTextEx(font, tmp, (Vector2){p1.x + (max_w.x - text_width.x / 2) + 20, p1.y - text_width.y / 2.0f}, font_size_small, spacing_small, RED);

        p1.x = min_x - 5.0f;
        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        // Grid line
        DrawLineV(p1, (Vector2){gui_state->ap_graph_config->max_x, p1.y}, LIGHTGRAY);

        // Tick line
        DrawLineV(p1, p2, RED);

        v += tick_ofsset;
    }

    // Draw vertical line
    {
        p1.x = min_x;
        p1.y = min_y + 5;

        p2.x = p1.x;
        p2.y = max_y;
        DrawLineV(p1, p2, RED);
    }

    // Draw horizontal line
    {
        p1.x = min_x;
        p1.y = min_y;

        p2.x = max_x;
        p2.y = p1.y;
        DrawLineV(p1, p2, RED);
    }

    struct action_potential *aps;

    // Draw the function
    for(int j = 0; j < n; j++) {

        aps = (struct action_potential *)gui_state->ap_graph_config->selected_aps[j].value;
        int c = arrlen(aps);
        int step = 1;

        if(c > 0) {
            Color line_color = colors[j % num_colors];
            for(int i = 0; i < c; i += step) {

                if(i + step < c && aps[i].t <= gui_config->time) {

                    p1.x = Remap(aps[i].t, 0.0f, gui_config->final_time, min_x, max_x);
                    p1.y = Remap(aps[i].v, gui_config->min_v, gui_config->max_v, min_y, max_y);

                    p2.x = Remap(aps[i + step].t, 0.0f, gui_config->final_time, min_x, max_x);
                    p2.y = Remap(aps[i + step].v, gui_config->min_v, gui_config->max_v, min_y, max_y);

                    if(aps[i + step].v > gui_config->max_v)
                        gui_config->max_v = aps[i + step].v;
                    if(aps[i + step].v < gui_config->min_v)
                        gui_config->min_v = aps[i + step].v;

                    DrawLineEx(p1, p2, 2.0f, line_color);
                }
            }
        }
    }

    // Draw AP coordinates over mouse cursor
    if(!gui_state->ap_graph_config->graph.drag && gui_state->ap_graph_config->selected_ap_point.x != FLT_MAX &&
       gui_state->ap_graph_config->selected_ap_point.y != FLT_MAX && gui_state->mouse_pos.x < max_x && gui_state->mouse_pos.x > min_x &&
       gui_state->mouse_pos.y < min_y && gui_state->mouse_pos.y > max_y) {

        char *tmp_point = "%.2lf, %.2lf";
        snprintf(tmp, TMP_SIZE, tmp_point, gui_state->ap_graph_config->selected_ap_point.x, gui_state->ap_graph_config->selected_ap_point.y);
        text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);
        DrawTextEx(font, tmp, (Vector2){gui_state->mouse_pos.x - text_width.x / 2, gui_state->mouse_pos.y - text_width.y}, font_size_small, spacing_small,
                   BLACK);
    }

    if(gui_state->ap_graph_config->selected_point_for_apd1.x != FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd1.y != FLT_MAX) {
        DrawCircleV(gui_state->ap_graph_config->selected_point_for_apd1, 4, RED);
    }

    if(gui_state->ap_graph_config->selected_point_for_apd2.x != FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd2.y != FLT_MAX) {
        DrawCircleV(gui_state->ap_graph_config->selected_point_for_apd2, 4, RED);
        DrawLineV(gui_state->ap_graph_config->selected_point_for_apd1, gui_state->ap_graph_config->selected_point_for_apd2, RED);

        float t1 = Remap(gui_state->ap_graph_config->selected_point_for_apd1.x, min_x, max_x, 0.0f, gui_config->final_time);
        float t2 = Remap(gui_state->ap_graph_config->selected_point_for_apd2.x, min_x, max_x, 0.0f, gui_config->final_time);

        char *tmp_point = "dt = %.4lf";
        snprintf(tmp, TMP_SIZE, tmp_point, fabsf(t2 - t1));
        text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

        float x = fminf(gui_state->ap_graph_config->selected_point_for_apd1.x, gui_state->ap_graph_config->selected_point_for_apd2.x);

        DrawTextEx(font, tmp, (Vector2){x + text_width.x / 2.0f, gui_state->ap_graph_config->selected_point_for_apd1.y - text_width.y}, font_size_small,
                   spacing_small, BLACK);
    }
}

void draw_scale(float min_v, float max_v, struct gui_state *gui_state, bool int_scale) {

    float scale_width = 20 * gui_state->ui_scale;

    float spacing_small = gui_state->font_spacing_small;
    float spacing_big = gui_state->font_spacing_big;

    uint32_t num_ticks;
    float tick_ofsset = 12.0f;

    if(!int_scale) {
        num_ticks = (uint32_t)((max_v - min_v) / tick_ofsset);

        if(num_ticks < MIN_SCALE_TICKS) {
            num_ticks = MIN_SCALE_TICKS;
            tick_ofsset = (max_v - min_v) / (float)num_ticks;
        } else if(num_ticks > MAX_SCALE_TICKS) {
            num_ticks = MAX_SCALE_TICKS;
            tick_ofsset = (max_v - min_v) / (float)num_ticks;
        }
    } else {
        num_ticks = 0;
        for(uint32_t i = (uint32_t)min_v; i < (uint32_t)max_v; i++) {
            num_ticks++;
        }
        tick_ofsset = (max_v - min_v) / (float)num_ticks;
    }

    char tmp[TMP_SIZE];

    float v = max_v;
    snprintf(tmp, TMP_SIZE, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(gui_state->font, tmp, gui_state->font_size_small, spacing_small);

    Vector2 p1, width;

    float scale_rec_height = 30.0f * gui_state->ui_scale;
    Color color;

    if(gui_state->scale.calc_bounds) {
        gui_state->scale.window.bounds.height += max_w.y;
        for(uint32_t t = 0; t <= num_ticks; t++) {
            gui_state->scale.window.bounds.height += scale_rec_height;
        }

        gui_state->scale.window.bounds.y -= gui_state->scale.window.bounds.height;

        gui_state->scale.calc_bounds = false;
    }

    Vector2 scale_bounds = (Vector2){(float)gui_state->scale.window.bounds.x, (float)gui_state->scale.window.bounds.y};

    if(!int_scale) {
        width = MeasureTextEx(gui_state->font, "Vm", gui_state->font_size_big, spacing_big);
        float diff = (float)scale_width - width.x;

        p1.x = scale_bounds.x + (diff / 2.0f);
        p1.y = scale_bounds.y - width.y;
        DrawTextEx(gui_state->font, "Vm", p1, gui_state->font_size_big, spacing_big, BLACK);
    }

    float initial_y = scale_bounds.y;

    for(uint32_t t = 0; t <= num_ticks; t++) {
        p1.x = scale_bounds.x - 10.0f * gui_state->ui_scale;
        p1.y = initial_y + (float)scale_rec_height / 2.0f;

        snprintf(tmp, TMP_SIZE, "%.2lf", v);
        width = MeasureTextEx(gui_state->font, tmp, gui_state->font_size_small, spacing_small);

        DrawTextEx(gui_state->font, tmp, (Vector2){p1.x - width.x, p1.y - width.y / 2.0f}, gui_state->font_size_small, spacing_small, BLACK);

        color = get_color((v - min_v) / (max_v - min_v), gui_state->scale_alpha, gui_state->current_scale);

        DrawRectangle((int)scale_bounds.x, (int)initial_y, (int)scale_width, (int)scale_rec_height, color);
        initial_y += scale_rec_height;
        v -= tick_ofsset;
    }
}

void draw_coordinates(struct gui_state *gui_state) {

    const float line_size = 1.0f;
    const float arrow_offset = 0.1f;
    static bool first_draw = true;

    if(first_draw) {
        gui_state->coordinates_cube = (Vector3){-(line_size / 2.0f) + 0.5f, -5.0f + 0.5f, -2.5f};
        first_draw = false;
    }

    Vector3 start_pos = (Vector3){gui_state->coordinates_cube.x - 0.5f, gui_state->coordinates_cube.y - 0.5f, gui_state->coordinates_cube.z - 0.5f};
    Vector3 end_pos = (Vector3){start_pos.x + line_size, start_pos.y, gui_state->coordinates_cube.z - 0.5f};

    DrawLine3D(start_pos, end_pos, RED);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y + arrow_offset, end_pos.z}, end_pos, RED);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, RED);

    gui_state->coordinates_label_x_position = GetWorldToScreen(end_pos, gui_state->camera);

    end_pos = (Vector3){start_pos.x, start_pos.y + line_size, end_pos.z};

    DrawLine3D(start_pos, end_pos, GREEN);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, GREEN);
    DrawLine3D((Vector3){end_pos.x + arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, GREEN);

    gui_state->coordinates_label_y_position = GetWorldToScreen((Vector3){end_pos.x, end_pos.y + 0.4f, end_pos.z}, gui_state->camera);

    end_pos = (Vector3){start_pos.x, start_pos.y, start_pos.z + line_size};

    DrawLine3D(start_pos, end_pos, DARKBLUE);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y, end_pos.z - arrow_offset}, end_pos, DARKBLUE);
    DrawLine3D((Vector3){end_pos.x + arrow_offset, end_pos.y, end_pos.z - arrow_offset}, end_pos, DARKBLUE);

    gui_state->coordinates_label_z_position = GetWorldToScreen(end_pos, gui_state->camera);
}

void reset_ui(struct gui_state *gui_state) {

    gui_state->help_box.window.bounds.x = 10;
    gui_state->help_box.window.bounds.y = 10;

    gui_state->ap_graph_config->graph.bounds.height = 300.0f * gui_state->ui_scale;
    gui_state->ap_graph_config->graph.bounds.width = 690.0f * gui_state->ui_scale;

    gui_state->ap_graph_config->graph.bounds.x = 10;
    gui_state->ap_graph_config->graph.bounds.y = (float)gui_state->current_window_height - gui_state->ap_graph_config->graph.bounds.height - 90;

    gui_state->search_window.bounds.width = 220;
    gui_state->search_window.bounds.height = 100;

    gui_state->mesh_info_box.window.bounds.x = (float)gui_state->current_window_width - gui_state->mesh_info_box.window.bounds.width - 10;
    gui_state->mesh_info_box.window.bounds.y = 10.0f;

    gui_state->end_info_box.window.bounds.x = gui_state->mesh_info_box.window.bounds.x - gui_state->mesh_info_box.window.bounds.width - 10;
    gui_state->end_info_box.window.bounds.y = gui_state->mesh_info_box.window.bounds.y;

    gui_state->scale.window.bounds.x = (float)gui_state->current_window_width - 30.0f * gui_state->ui_scale;
    gui_state->scale.window.bounds.y = (float)gui_state->current_window_height / 1.5f;

    gui_state->scale.window.bounds.width = 20;
    gui_state->scale.window.bounds.height = 0;
    gui_state->scale.calc_bounds = true;

    gui_state->controls_window.bounds.x = (float)gui_state->current_window_width / 2.0f - gui_state->controls_window.bounds.width;
    gui_state->controls_window.bounds.y = 10;
}

void reset(struct gui_shared_info *gui_config, struct gui_state *gui_state, bool full_reset) {

    if(!gui_config->paused) {
        return;
    }

    gui_config->grid_info.loaded = false;

    gui_state->voxel_alpha = 255;

    size_t n = hmlen(gui_state->ap_graph_config->selected_aps);

    for(size_t i = 0; i < n; i++) {
        arrsetlen(gui_state->ap_graph_config->selected_aps[i].value, 0);
    }

    gui_config->restart = true;

    if(gui_config->draw_type == DRAW_SIMULATION) {
        gui_config->grid_info.alg_grid = NULL;
    } else {
        free_vtk_unstructured_grid(gui_config->grid_info.vtk_grid);
        gui_config->grid_info.vtk_grid = NULL;
    }

    gui_state->ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    omp_unset_nest_lock(&gui_config->sleep_lock);
    gui_state->current_selected_volume.position_draw = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->current_mouse_over_volume.position_draw = (Vector3){-1, -1, -1};
    gui_state->ap_graph_config->selected_point_for_apd1 = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd2 = (Vector2){FLT_MAX, FLT_MAX};

    if(full_reset) {

        for(size_t i = 0; i < n; i++) {
            arrfree(gui_state->ap_graph_config->selected_aps[i].value);
        }

        hmfree(gui_state->ap_graph_config->selected_aps);
        gui_state->ap_graph_config->selected_aps = NULL;
        hmdefault(gui_state->ap_graph_config->selected_aps, NULL);

        hmfree(gui_state->current_selected_volumes);
        gui_state->current_selected_volumes = NULL;

        set_camera_params(&(gui_state->camera), false);

        gui_state->scale_alpha = 255;

        reset_ui(gui_state);

        gui_state->show_coordinates = true;
        gui_state->slicing_mesh = false;
        gui_state->current_mode = VISUALIZING;

        gui_state->plane_roll = 0.0f;
        gui_state->plane_pitch = 0.0f;
        gui_state->plane_tx = 0.0f;
        gui_state->plane_ty = 0.0f;
        gui_state->plane_tz = 0.0f;
    }
}

void draw_control_window(struct gui_state *gui_state, struct gui_shared_info *gui_config) {

    static bool spinner_edit = false;
    gui_state->controls_window.show = !GuiWindowBox(gui_state->controls_window.bounds, "Controls");

    Rectangle button_pos = (Rectangle){0, 0, 32, 32};

    button_pos.x = gui_state->controls_window.bounds.x + 2.0f;
    button_pos.y = gui_state->controls_window.bounds.y + WINDOW_STATUSBAR_HEIGHT + 3.0f;

    bool update_main = false;

    DISABLE_IF_NOT_PAUSED

    if(GuiButton(button_pos, "#76#")) {
        reset(gui_config, gui_state, false);
    }
    ENABLE;

    DISABLE_IF_NOT_PAUSED_OR_NOT_IN_DRAW
    DISABLE_IF_IN_DRAW_AND_SINGLE_FILE

    button_pos.x += button_pos.width + 4.0f;

    if(GuiButton(button_pos, "#129#")) {
        gui_config->current_file_index = 0;
        update_main = true;
    }

    // DISABLE_IF_IN_DRAW_AND_SINGLE_FILE
    //  return button
    {
        button_pos.x += button_pos.width + 4.0f;

        if(GuiButton(button_pos, "#114#")) {
            if(gui_config->paused) {
                gui_config->current_file_index -= 1;
                update_main = true;
            }
        }
    }

    ENABLE;

    DISABLE_IF_IN_DRAW_AND_SINGLE_FILE

    // Play or pause button
    {
        button_pos.x += button_pos.width + 4.0f;

        if(gui_config->paused) {
            gui_config->paused = !GuiButton(button_pos, "#131#");
        } else {
            gui_config->paused = GuiButton(button_pos, "#132#");
        }
    }

    DISABLE_IF_NOT_PAUSED_OR_NOT_IN_DRAW
    // advance button
    {
        button_pos.x += button_pos.width + 4.0f;

        if(GuiButton(button_pos, "#115#")) {
            if(gui_config->paused) {
                gui_config->current_file_index += 1;
                update_main = true;
            }
        }
    }

    button_pos.x += button_pos.width + 4.0f;
    if(GuiButton(button_pos, "#134#")) {
        update_main = true;
        gui_config->current_file_index = (float)gui_config->final_file_index;
    }

    button_pos.x += button_pos.width + 4.0f;

    ENABLE;
    DISABLE_IF_NOT_PAUSED_OR_NOT_IN_DRAW
    DISABLE_IF_IN_DRAW_AND_SINGLE_FILE
    // Calc bounds button
    if(GuiButton(button_pos, "#94#")) {
        gui_config->calc_bounds = true;
        omp_unset_nest_lock(&gui_config->sleep_lock);
    }

    button_pos.x += button_pos.width + 4.0f;
    button_pos.width = 96;

    ENABLE;
    DISABLE_IF_NOT_PAUSED_OR_NOT_IN_DRAW
    DISABLE_IF_IN_DRAW_AND_SINGLE_FILE

    int old_index = (int)gui_config->current_file_index;

    if(NOT_IN_DRAW) {
        GuiSpinner(button_pos, NULL, &gui_config->time, 0, gui_config->final_time, false, spinner_edit);
    } else if(GuiSpinner(button_pos, NULL, &gui_config->current_file_index, 0, (float)gui_config->final_file_index, true, spinner_edit)) {
        spinner_edit = !spinner_edit;
    }

    if(old_index != (int)gui_config->current_file_index) {
        update_main = true;
    }
    gui_state->handle_keyboard_input = !spinner_edit;

    ENABLE;

    if(update_main && gui_config->draw_type == DRAW_FILE) {

        CHECK_FILE_INDEX(gui_config)

        if(gui_config->current_file_index == 0) {
            size_t n = hmlen(gui_state->ap_graph_config->selected_aps);
            for(size_t i = 0; i < n; i++) {
                arrsetlen(gui_state->ap_graph_config->selected_aps[i].value, 0);
            }
        }

        omp_unset_nest_lock(&gui_config->sleep_lock);
    }
}
