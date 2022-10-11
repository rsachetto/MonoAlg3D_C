SOURCE_FILES="gui.c gui_colors.c gui_window_helpers.c gui_mesh_helpers.c gui_draw.c raylib_ext.c"
HEADER_FILES="gui.h gui_colors.h gui_window_helpers.h gui_mesh_helpers.h gui_draw.h raylib_ext.h ../3dparty/raylib/src/extras/raygui.h"
COMPILE_STATIC_LIB "gui" "$SOURCE_FILES" "$HEADER_FILES"
