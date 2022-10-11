#ifndef __GUI_WINDOW_HELPERS_H
#define __GUI_WINDOW_HELPERS_H

#include "gui.h"
#include "../3dparty/raylib/src/extras/raygui.h"
#include "../3dparty/raylib/src/extras/gui_textbox_extended.h"

void check_colisions_for_move(struct gui_state *gui_state);
void maybe_move_or_drag(struct gui_state *gui_state);
void draw_text_window(struct gui_text_window *box, struct gui_state *gui_state, float text_offset);
bool draw_search_window(struct gui_state *gui_state);

#endif /* __GUI_WINDOW_HELPERS_H */
