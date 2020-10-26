//
// Created by sachetto on 11/11/17.
//

#ifndef MONOALG3D_GUI_H
#define MONOALG3D_GUI_H

#include "../alg/grid/grid.h"
#include "../alg/cell/cell.h"
#include "../vtk_utils/vtk_unstructured_grid.h"
#include "../3dparty/raylib/src/raymath.h"
#include "../config/config_parser.h"

#include <omp.h>

#define MIN_VERTICAL_TICKS 4
#define MAX_VERTICAL_TICKS 20

#define MIN_HORIZONTAL_TICKS 4
#define MAX_HORIZONTAL_TICKS 20

#define SIZEOF(A) (sizeof(A)/sizeof(A[0]))
#define WIDER_TEXT "------------------------------------------------------"
#define DOUBLE_CLICK_DELAY 0.8 //seconds

struct action_potential {
    real_cpu v;
    real_cpu t;
};

enum draw_simulation_or_file {
    DRAW_SIMULATION,
    DRAW_FILE,
};

struct vector3_voidp_hash_entry {
    Vector3 key;
    void *value;
};


struct ap_graph_config {

    Rectangle graph;
    float max_x;
    float min_y;
    float max_y;
    float min_x;

    Vector2 selected_ap_point;
    Vector2 selected_point_for_apd1;
    Vector2 selected_point_for_apd2;
    struct vector3_voidp_hash_entry *selected_aps;
    Rectangle drag_graph_button_position;
    Rectangle move_graph_button_position;
    bool drag_ap_graph;
    bool move_ap_graph;
    bool draw_selected_ap_text;

};

typedef struct action_potential * action_potential_array;

struct voxel {
	Vector3 position_draw;
	Vector3 position_mesh;
	Vector3 size;
	real_cpu v;
	uint32_t matrix_position;
};

struct gui_config {

    real_cpu max_v;
    real_cpu min_v;
    bool simulating;
    bool paused;
    bool exit;
    bool adaptive;
    bool restart;
    real_cpu time;
    real_cpu final_time;
    real_cpu dt;
    int step;

    char *config_name;
    char *input;

    long solver_time;
    long ode_total_time;
    long cg_total_time;
    long total_mat_time;
    long total_ref_time;
    long total_deref_time;
    long total_write_time;
    long total_config_time;
    long total_cg_it;

    int advance_or_return;
    int draw_type;
    //If we are compiling this file, openmp is available.
    omp_lock_t draw_lock;
    omp_lock_t sleep_lock;

    struct grid_info {
        union {
            struct vtk_unstructured_grid *vtk_grid;
            struct grid *alg_grid;
        };
        char *file_name;
        bool loaded;
    } grid_info;

    char *error_message;
    bool int_scale;
};

struct gui_state {

    bool handle_keyboard_input;
    bool one_selected;
    bool show_ap;
    bool c_pressed;
    bool draw_grid_lines;
    bool draw_grid_only ;

    Rectangle help_box;
    bool show_help_box;
    bool move_help_box;

    Rectangle mesh_info_box;
    bool show_mesh_info_box;
    bool move_info_box;

    bool show_selection_box;

    Rectangle scale_bounds;
    bool show_scale;
    bool move_scale;

    Rectangle end_info_box;
    bool show_end_info_box;
    bool move_end_info_box;

    Ray ray;
	Ray ray_mouse_over;

    int current_scale;
    int scale_alpha;
    int voxel_alpha;

    struct voxel current_selected_volume;
    struct voxel current_mouse_over_volume;

    //Font font;
    float font_size_small;
    float font_size_big;

    double mouse_timer;
    double selected_time;

    Vector2 mouse_pos;
    Vector2 sub_window_pos;
    bool move_sub_window;

    float box_width;
    float box_height;

    struct ap_graph_config *ap_graph_config;

    Vector3 coordinates_cube;
    Vector3 coordinates_cube_size;

    bool show_coordinates;
    bool move_coordinates;

    Vector2 coordinates_label_x_position;
    Vector2 coordinates_label_y_position;
    Vector2 coordinates_label_z_position;

    Font font;

    Camera3D camera;

    bool double_clicked;

};

struct mesh_info {
    bool center_calculated;
    Vector3 max_size;
    Vector3 min_size;
};

void init_and_open_gui_window();
void gui_set_alg_grid(struct grid *the_grid);

void gui_set_vtk_grid(struct vtk_unstructured_grid *the_grid);
struct vtk_unstructured_grid * gui_get_vtk_grid();

void gui_set_simulating(bool state);
void gui_set_paused(bool paused);
bool gui_get_paused();
void gui_set_grid_loaded(bool loaded);
void gui_lock_sleep_lock();
void gui_unlock_sleep_lock();
void gui_lock_draw_lock();
void gui_unlock_draw_lock();
char * gui_get_error_message();
void gui_free_error_message();
void gui_set_error_message(char *msg);
real_cpu gui_get_dt();
void gui_set_dt(real_cpu dt);

int gui_get_step();
void gui_set_step(int s);

char * gui_get_filename();
void gui_set_filename(char *filename);

real_cpu gui_get_final_time();
void gui_set_final_time(real_cpu time);

int gui_get_advance_or_return();
void gui_set_advance_or_return(int adv);

char * gui_get_input();
void gui_set_input(char *input);

void gui_end_simulation(long res_time, long ode_total_time, long cg_total_time, long total_mat_time, long total_ref_time, long total_deref_time, long total_write_time, long total_config_time, long total_cg_it);
void init_gui_config_for_simulation(struct user_options *options);
void read_and_render_activation_map(char *input_file, char *error);
void init_gui_config_for_visualization(struct visualization_options *options, bool only_restart);

bool gui_get_restart();
bool gui_get_exit();
void gui_set_time(double time);

#endif //MONOALG3D_DRAW_H
