#include "gui_window_helpers.h"
#include <float.h>

static inline void move_rect(Vector2 new_pos, Rectangle *rect) {

    float new_x = new_pos.x;
    float new_y = new_pos.y;

    if(new_x > 1 && new_x + rect->width < (float)GetScreenWidth()) {
        rect->x = new_x;
    }
    if(new_y > 1 && new_y + rect->height < (float)GetScreenHeight()) {
        rect->y = new_y;
    }
}

static inline void drag_box(Vector2 mouse_pos, Rectangle *box) {
    float new_x = mouse_pos.x + box->width / 2;
    float new_y = mouse_pos.y - box->height / 2;
    move_rect((Vector2){new_x, new_y}, box);
}

#define COLIDE_STATUS_BAR_WITH_OFFSET(mouse_pos, bounds, y_off)                                                                                                \
    CheckCollisionPointRec(mouse_pos, (Rectangle){bounds.x, bounds.y - y_off, bounds.width - 18, WINDOW_STATUSBAR_HEIGHT})
#define COLIDE_STATUS_BAR(mouse_pos, bounds) COLIDE_STATUS_BAR_WITH_OFFSET(mouse_pos, bounds, 0)

void check_colisions_for_move(struct gui_state *gui_state) {

    if(COLIDE_STATUS_BAR(gui_state->mouse_pos, gui_state->search_window.bounds)) {
        gui_state->search_window.move = true;
    } else if(COLIDE_STATUS_BAR(gui_state->mouse_pos, gui_state->help_box.window.bounds)) {
        gui_state->help_box.window.move = true;
    } else if(COLIDE_STATUS_BAR(gui_state->mouse_pos, gui_state->slice_help_box.window.bounds)) {
        gui_state->slice_help_box.window.move = true;
    } else if(COLIDE_STATUS_BAR(gui_state->mouse_pos, gui_state->mesh_info_box.window.bounds)) {
        gui_state->mesh_info_box.window.move = true;
    } else if(COLIDE_STATUS_BAR(gui_state->mouse_pos, gui_state->end_info_box.window.bounds)) {
        gui_state->end_info_box.window.move = true;
    } else if(COLIDE_STATUS_BAR(gui_state->mouse_pos, gui_state->controls_window.bounds)) {
        gui_state->controls_window.move = true;
    } else if(COLIDE_STATUS_BAR_WITH_OFFSET(gui_state->mouse_pos, gui_state->ap_graph_config->graph.bounds, WINDOW_STATUSBAR_HEIGHT)) {
        gui_state->ap_graph_config->graph.move = true;
    } else if(CheckCollisionPointRec(gui_state->mouse_pos, (Rectangle){gui_state->scale.window.bounds.x, gui_state->scale.window.bounds.y,
                                                                       gui_state->scale.window.bounds.width, gui_state->scale.window.bounds.height})) {
        gui_state->scale.window.move = true;
    } else if(CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->drag_graph_button)) {
        gui_state->ap_graph_config->graph.drag = true;
    }
}

#define PERFORM_MOVE(mouse_pos, win)                                                                                                                           \
    do {                                                                                                                                                       \
        move_rect((Vector2){(mouse_pos.x) - (win.bounds.width - 18.0f) / 2.0f, (mouse_pos.y) + WINDOW_STATUSBAR_HEIGHT / 2.0f}, &(win.bounds));                \
        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON))                                                                                                           \
            win.move = false;                                                                                                                                  \
    } while(0)

void maybe_move_or_drag(struct gui_state *gui_state) {

    if(gui_state->search_window.move) {
        gui_state->search_window.bounds.x = (gui_state->mouse_pos.x) - (gui_state->search_window.bounds.width - 18.0f) / 2.0f;
        gui_state->search_window.bounds.y = (gui_state->mouse_pos.y) - WINDOW_STATUSBAR_HEIGHT / 2.0f;

        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) {
            gui_state->search_window.move = false;
        }
    } else if(gui_state->ap_graph_config->graph.drag) {

        float new_heigth = gui_state->mouse_pos.y - gui_state->ap_graph_config->graph.bounds.y;

        if(new_heigth > 100) {
            gui_state->ap_graph_config->graph.bounds.height = new_heigth;
        }

        float new_width = gui_state->mouse_pos.x - gui_state->ap_graph_config->graph.bounds.x;

        if(new_width > 200) {
            gui_state->ap_graph_config->graph.bounds.width = new_width;
        }

        gui_state->ap_graph_config->selected_point_for_apd1.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd1.y = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;

        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) {
            gui_state->ap_graph_config->graph.drag = false;
        }
    } else if(gui_state->ap_graph_config->graph.move) {

        PERFORM_MOVE(gui_state->mouse_pos, gui_state->ap_graph_config->graph);

        gui_state->ap_graph_config->selected_point_for_apd1.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd1.y = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;

    } else if(gui_state->help_box.window.move) {
        PERFORM_MOVE(gui_state->mouse_pos, gui_state->help_box.window);
    } else if(gui_state->slice_help_box.window.move) {
        PERFORM_MOVE(gui_state->mouse_pos, gui_state->slice_help_box.window);
    } else if(gui_state->mesh_info_box.window.move) {
        PERFORM_MOVE(gui_state->mouse_pos, gui_state->mesh_info_box.window);
    } else if(gui_state->end_info_box.window.move) {
        PERFORM_MOVE(gui_state->mouse_pos, gui_state->end_info_box.window);
    } else if(gui_state->controls_window.move) {
        PERFORM_MOVE(gui_state->mouse_pos, gui_state->controls_window);
    }
    else if(gui_state->scale.window.move) {
        drag_box(gui_state->mouse_pos, &gui_state->scale.window.bounds);
        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON))
            gui_state->scale.window.move = false;
    }
}

void draw_text_window(struct gui_text_window *box, struct gui_state *gui_state, float text_offset) {

    float text_x = box->window.bounds.x + 20;
    float text_y = box->window.bounds.y + 10 + WINDOW_STATUSBAR_HEIGHT;

    box->window.show = !GuiWindowBox(box->window.bounds, box->title);

    int num_lines = box->num_lines;

    for(int i = 0; i < num_lines; i++) {
        DrawTextEx(gui_state->font, box->lines[i], (Vector2){text_x, text_y}, gui_state->font_size_small, gui_state->font_spacing_small, BLACK);
        text_y += text_offset;
    }
}

bool draw_search_window(struct gui_state *gui_state) {

#define CENTER_X "Center X"
#define CENTER_Y "Center Y"
#define CENTER_Z "Center Z"

    Vector2 text_box_size = MeasureTextEx(gui_state->font, CENTER_X, gui_state->font_size_small, gui_state->font_spacing_small);

    float text_box_y_dist = text_box_size.y * 3.0f;
    float label_box_y_dist = 30;
    float x_off = 10;

    static char center_x_text[128] = {0};
    static char center_y_text[128] = {0};
    static char center_z_text[128] = {0};

    float pos_x = gui_state->search_window.bounds.x;
    float pos_y = gui_state->search_window.bounds.y;

    float box_pos = pos_x + x_off;
    gui_state->search_window.bounds.width = text_box_size.x * 3.8f;
    gui_state->search_window.bounds.height = (text_box_size.y + text_box_y_dist) * 1.6f;

    bool window_closed = GuiWindowBox(gui_state->search_window.bounds, "Enter the center of the cell");

    DrawTextEx(gui_state->font, CENTER_X, (Vector2){box_pos, pos_y + label_box_y_dist}, gui_state->font_size_small, gui_state->font_spacing_small, BLACK);

    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_size.x, text_box_size.y}, center_x_text, SIZEOF(center_x_text) - 1, true);

    box_pos = pos_x + text_box_size.x + 2 * x_off;
    DrawTextEx(gui_state->font, CENTER_Y, (Vector2){box_pos, pos_y + label_box_y_dist}, gui_state->font_size_small, gui_state->font_spacing_small, BLACK);
    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_size.x, text_box_size.y}, center_y_text, SIZEOF(center_y_text) - 1, true);

    box_pos = pos_x + 2 * text_box_size.x + 3 * x_off;
    DrawTextEx(gui_state->font, CENTER_Z, (Vector2){box_pos, pos_y + label_box_y_dist}, gui_state->font_size_small, gui_state->font_spacing_small, BLACK);
    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_size.x, text_box_size.y}, center_z_text, SIZEOF(center_z_text) - 1, true);

    bool btn_ok_clicked =
        GuiButton((Rectangle){pos_x + text_box_size.x + 2 * x_off, pos_y + (text_box_size.y + text_box_y_dist) * 1.2f, text_box_size.x, text_box_size.y}, "OK");

    if(btn_ok_clicked) {
        gui_state->found_volume.position_mesh.x = strtof(center_x_text, NULL);
        gui_state->found_volume.position_mesh.y = strtof(center_y_text, NULL);
        gui_state->found_volume.position_mesh.z = strtof(center_z_text, NULL);
    }

    return window_closed || btn_ok_clicked;
}
