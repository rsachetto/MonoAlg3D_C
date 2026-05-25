//
// Created by sachetto on 11/11/17.
//

#ifndef MONOALG3D_GUI_H
#define MONOALG3D_GUI_H

#include <omp.h>
#include <sys/types.h>

#include "../3dparty/raylib/src/raylib.h"
#include "../3dparty/raylib/src/extras/rlights.h"
#include "../3dparty/raylib/src/raymath.h"
#include "../alg/cell/cell.h"
#include "../alg/grid/grid.h"
#include "../config/config_parser.h"
#include "../vtk_utils/vtk_unstructured_grid.h"

#include "gui_colors.h"

#define MIN_VERTICAL_TICKS 4
#define MAX_VERTICAL_TICKS 20

#define MIN_HORIZONTAL_TICKS 4
#define MAX_HORIZONTAL_TICKS 20

#define SIZEOF(A) (sizeof(A) / sizeof(A[0]))
#define WIDER_TEXT " - Alt + R to restart simulation and the box positions"
#define DOUBLE_CLICK_DELAY 0.7 // seconds

#define TMP_SIZE 256

#define V3_SAME(v)                                                                                                                                             \
    (Vector3) {                                                                                                                                                \
        v, v, v                                                                                                                                                \
    }

#define WINDOW_STATUSBAR_HEIGHT 22

#define CHECK_FILE_INDEX(gui_config)                                                                                                                           \
    if((gui_config)->current_file_index < 0)                                                                                                                   \
        (gui_config)->current_file_index++;                                                                                                                    \
    else if((gui_config)->current_file_index > (gui_config)->final_file_index)                                                                                 \
        (gui_config)->current_file_index = (gui_config)->final_file_index;

#define NOT_PAUSED !gui_config->paused
#define NOT_IN_DRAW gui_config->draw_type != DRAW_FILE
#define IN_DRAW gui_config->draw_type == DRAW_FILE

#define DISABLE_IF_NOT_PAUSED                                                                                                                                  \
    if(NOT_PAUSED) {                                                                                                                                           \
        GuiDisable();                                                                                                                                          \
    }

#define DISABLE_IF_NOT_PAUSED_OR_NOT_IN_DRAW                                                                                                                   \
    if(NOT_PAUSED || NOT_IN_DRAW) {                                                                                                                            \
        GuiDisable();                                                                                                                                          \
    }

#define ENABLE GuiEnable()

#define DISABLE_IF_IN_DRAW_AND_SINGLE_FILE                                                                                                                     \
    if(gui_config->draw_type == DRAW_FILE && gui_config->final_file_index == 0) {                                                                              \
        GuiDisable();                                                                                                                                          \
    }

struct action_potential {
    float v;
    float t;
};

enum draw_simulation_or_file {
    DRAW_SIMULATION,
    DRAW_FILE,
};

enum mode {
    VISUALIZING,
    SLICING,
    EDITING,
};

struct vector3_voidp_hash_entry {
    Vector3 key;
    void *value;
};

typedef struct action_potential *action_potential_array;

struct voxel {
    Vector3 position_draw;
    Vector3 position_mesh;
    Vector3 size;
    uint32_t draw_index;
    uint32_t matrix_position;
    float v;
};

struct vector3_voxel_entry {
    Vector3 key;
    struct voxel value;
};

// This struct is shared with the main thread
struct gui_shared_info {
    float max_v;
    float min_v;
    bool simulating;
    bool paused;
    bool exit;
    bool adaptive;
    bool restart;
    bool enable_slice;
    float time;
    float final_time;
    float current_file_index;
    int final_file_index;
    float dt;
    int step;
    size_t progress;
    size_t file_size;

    char *config_name;
    char *input;

    uint64_t solver_time;
    uint64_t ode_total_time;
    uint64_t cg_total_time;
    uint64_t total_mat_time;
    uint64_t total_ref_time;
    uint64_t total_deref_time;
    uint64_t total_write_time;
    uint64_t total_cg_it;

    int draw_type;

    // If we are compiling this file, openmp is available.
    omp_nest_lock_t draw_lock;
    omp_nest_lock_t sleep_lock;

    struct grid_info {
        union {
            struct vtk_unstructured_grid *vtk_grid;
            struct grid *alg_grid;
        };
        char *file_name;
        bool loaded;
    } grid_info;

    char *message;
    bool int_scale;
    float ui_scale;

    struct simulation_files *simulation_files;

    bool calc_bounds;
};

struct window_commom {
    Rectangle bounds;
    bool show;
    bool move;
    bool drag;
};

struct gui_text_window {
    struct window_commom window;
    char **lines;
    char *title;
    int num_lines;
};

struct gui_scale {
    struct window_commom window;
    bool calc_bounds;
};

struct ap_graph_config {

    struct window_commom graph;
    float max_x;
    float min_y;
    float max_y;
    float min_x;

    Vector2 selected_ap_point;
    Vector2 selected_point_for_apd1;
    Vector2 selected_point_for_apd2;
    struct vector3_voidp_hash_entry *selected_aps;
    Rectangle drag_graph_button;
};

struct gui_state {

    int current_window_width;
    int current_window_height;

    int current_data_index;
    int max_data_index;

    float ui_scale;

    float font_spacing_big;
    float font_spacing_small;

    bool handle_keyboard_input;
    bool c_pressed;
    bool draw_grid_lines;
    bool draw_grid_only;

    struct gui_text_window help_box;
    struct gui_text_window slice_help_box;
    struct gui_text_window edit_help_box;
    struct gui_text_window mesh_info_box;
    struct gui_text_window end_info_box;
    struct window_commom search_window;
    struct window_commom controls_window;

    struct gui_scale scale;

    Ray ray;
    float ray_hit_distance;
    Ray ray_mouse_over;
    float ray_mouse_over_hit_distance;

    int current_scale;
    int scale_alpha;
    int voxel_alpha;

    struct vector3_voxel_entry *current_selected_volumes;
    struct voxel current_selected_volume;
    struct voxel found_volume;
    struct voxel current_mouse_over_volume;

    float font_size_small;
    float font_size_big;

    double mouse_timer;

    Vector2 mouse_pos;

    struct ap_graph_config *ap_graph_config;

    Vector3 coordinates_cube;
    bool show_coordinates;

    Vector2 coordinates_label_x_position;
    Vector2 coordinates_label_y_position;
    Vector2 coordinates_label_z_position;

    Font font;

    Camera3D camera;

    Light light;

    bool double_clicked;

    Vector3 plane_normal;
    Vector3 plane_point;

    enum mode current_mode;

    bool slicing_mesh;
    bool recalculating_visibility;

    bool plane_loaded;
    bool visibility_recalculated;
    ui8_array old_cell_visibility;
    bool *exclude_from_mesh;

    bool ctrl_pressed;

    float plane_roll;
    float plane_pitch;
    float plane_tx;
    float plane_ty;
    float plane_tz;

    float mesh_scale_factor;
    Vector3 mesh_offset;
    u_int32_t last_n_active;

    bool get_cell_property;
    bool paste_cell_property;
    float copied_property_value;
};

struct mesh_info {
    bool center_calculated;
    Vector3 max_size;
    Vector3 min_size;
};

#ifdef __cplusplus
extern "C" {
#endif

void init_and_open_gui_window(struct gui_shared_info *gui_config);

#ifdef __cplusplus
}
#endif

#endif // MONOALG3D_GUI_H
