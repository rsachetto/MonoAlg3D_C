
#ifndef MONOALG3D_GUI_DRAW_H
#define MONOALG3D_GUI_DRAW_H

#include "gui.h"
#include "../3dparty/raylib/src/extras/raygui.h"

#define MIN_SCALE_TICKS 5
#define MAX_SCALE_TICKS 12

void set_camera_params(Camera3D *camera, bool set_mode);
void draw_ap_graph(struct gui_state *gui_state, struct gui_shared_info *gui_config);
void draw_scale(float min_v, float max_v, struct gui_state *gui_state, bool int_scale);
void draw_coordinates(struct gui_state *gui_state);
void reset_ui(struct gui_state *gui_state);
void reset(struct gui_shared_info *gui_config, struct gui_state *gui_state, bool full_reset);
void draw_control_window(struct gui_state *gui_state, struct gui_shared_info *gui_config);

#endif //MONOALG3D_GUI_DRAW_H
