//
// Created by sachetto on 11/11/17.
//
#include <float.h>
#include <time.h>
#include <string.h>

#include "../3dparty/tinyfiledialogs/tinyfiledialogs.h"
#include "../logger/logger.h"
#include "../3dparty/stb_ds.h"
#include "../utils/file_utils.h"
#include "color_maps.h"
#include "gui.h"

//RAYLIB//////
#include "../3dparty/raylib/ricons.h"
#include "../3dparty/raylib/src/camera.h"

#define RAYGUI_IMPLEMENTATION
#include "../3dparty/raylib/src/raygui.h"
#undef RAYGUI_IMPLEMENTATION

#define GUI_TEXTBOX_EXTENDED_IMPLEMENTATION
#include "../3dparty/raylib/src/external/glad.h"
#include "../3dparty/raylib/src/gui_textbox_extended.h"
#include "../3dparty/raylib/src/rlgl.h"
/////////

static int current_window_height = 0;
static int current_window_width = 0;

int info_box_lines;
const int end_info_box_lines = 10;
const int mesh_info_box_lines = 9;

void set_camera_params(Camera3D *camera) {
    camera->position = (Vector3){ 0.1f, 0.1f, 20.f};  // Camera position
    camera->target = (Vector3){ 0.f, 0.f, 0.f};
    camera->up       = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera->fovy     = 45.0f;                                // Camera field-of-view Y
    camera->type     = CAMERA_PERSPECTIVE;                  // Camera mode type
    SetCameraMode(*camera, CAMERA_FREE); // Set a free camera mode
}

static struct gui_state * new_gui_state_with_font_sizes(int font_size_small, int font_size_big) {

    struct gui_state *gui_state = (struct gui_state*) calloc(1, sizeof(struct gui_state));

    set_camera_params(&(gui_state->camera));

    gui_state->font = GetFontDefault();

    gui_state->handle_keyboard_input = true;
    gui_state->font_size_small = font_size_small;
    gui_state->font_size_big = font_size_big;

    gui_state->show_ap = true;
    gui_state->show_scale = true;

    gui_state->show_help_box = true;
    gui_state->help_box.x = 10;
    gui_state->help_box.y = 10;

    gui_state->show_end_info_box = true;
    gui_state->show_mesh_info_box = true;
    gui_state->ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->current_selected_volume = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->current_mouse_over_volume.position_draw = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->mouse_pos = (Vector2){0, 0};
    gui_state->current_scale = 0;
    gui_state->voxel_alpha = 255;
    gui_state->scale_alpha = 255;
    gui_state->mouse_timer = -1;
    gui_state->selected_time = 0.0;
    gui_state->move_sub_window = false;

    gui_state->ap_graph_config = (struct ap_graph_config*) malloc(sizeof(struct ap_graph_config));
    gui_state->ap_graph_config->selected_ap_point = (Vector2) {FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd1 = (Vector2) {FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd2 = (Vector2) {FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_aps = NULL;
    gui_state->ap_graph_config->drag_ap_graph = false;
    gui_state->ap_graph_config->move_ap_graph = false;

    gui_state->box_width = 220;
    gui_state->box_height = 100;

    hmdefault(gui_state->ap_graph_config->selected_aps, NULL);

    gui_state->ap_graph_config->graph.height = (float)450;
    gui_state->ap_graph_config->graph.width = (float)900;

    gui_state->ap_graph_config->graph.x = 100;
    gui_state->ap_graph_config->graph.y = (float)current_window_height - gui_state->ap_graph_config->graph.height - 70;
    gui_state->scale_bounds.x = (float)current_window_width - 30.0f;
    gui_state->scale_bounds.y = (float)current_window_height/1.5f;

    gui_state->scale_bounds.width = 20;
    gui_state->scale_bounds.height = 0;

    gui_state->show_coordinates = true;

    gui_state->coordinates_cube_size = (Vector3){1.2,1.2,1.2};

    gui_state->double_clicked = false;

    return gui_state;
}

struct mesh_info *new_mesh_info() {
    struct mesh_info *mesh_info = (struct mesh_info *) malloc(sizeof(struct mesh_info));
    mesh_info->center_calculated = false;
    return mesh_info;
}

static inline float normalize(float r_min, float r_max, float t_min, float t_max, float m) {
    return ((m - r_min) / (r_max-r_min))*(t_max - t_min) + t_min;
}

static inline Color get_color(real_cpu value, int alpha, int current_scale) {

    int idx1;        // |-- Our desired color will be between these two indexes in "color".
    int idx2;        // |
    real_cpu fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

    if(value <= 0)      {  idx1 = idx2 = 0;            }    // accounts for an input <=0
    else if(value >= 1)  {  idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
    else
    {
        value = value * (NUM_COLORS-1);        // Will multiply value by NUM_COLORS.
        idx1  = (int)floor(value);                  // Our desired color will be after this index.
        idx2  = idx1+1;                        // ... and before this index (inclusive).
        fractBetween = value - (real_cpu)idx1;    // Distance between the two indexes (0-1).
    }

    unsigned char red;
    unsigned char green;
    unsigned char blue;

    red  = (unsigned char) (((color_scales[current_scale][idx2][0] - color_scales[current_scale][idx1][0])*fractBetween + color_scales[current_scale][idx1][0]) * 255);
    green = (unsigned char) (((color_scales[current_scale][idx2][1] - color_scales[current_scale][idx1][1])*fractBetween + color_scales[current_scale][idx1][1]) * 255);
    blue  = (unsigned char) (((color_scales[current_scale][idx2][2] - color_scales[current_scale][idx1][2])*fractBetween + color_scales[current_scale][idx1][2]) * 255);

    Color result;
    result.r = red;
    result.g = green;
    result.b = blue;
    result.a = alpha;

    return result;
}

static Vector3 find_mesh_center(struct mesh_info *mesh_info) {

    struct grid *grid_to_draw = gui_config.grid_info.alg_grid;

    uint32_t n_active = grid_to_draw->num_active_cells;
    struct cell_node **ac = grid_to_draw->active_cells;
    struct cell_node *grid_cell;

    float mesh_max_x, mesh_max_y, mesh_max_z;
    float mesh_min_x, mesh_min_y, mesh_min_z;

    mesh_max_x = FLT_MIN;
    mesh_max_y = FLT_MIN;
    mesh_max_z = FLT_MIN;

    mesh_min_x = FLT_MAX;
    mesh_min_y = FLT_MAX;
    mesh_min_z = FLT_MAX;

    float mesh_max_dx, mesh_max_dy, mesh_max_dz;
    float mesh_min_dx, mesh_min_dy, mesh_min_dz;

    mesh_max_dx = FLT_MIN;
    mesh_max_dy = FLT_MIN;
    mesh_max_dz = FLT_MIN;

    mesh_min_dx = FLT_MAX;
    mesh_min_dy = FLT_MAX;
    mesh_min_dz = FLT_MAX;

    Vector3 result = (Vector3){0.0f, 0.0f, 0.0f};

    if (ac) {
        for (uint32_t i = 0; i < n_active; i++) {
            grid_cell = ac[i];
            if(grid_cell->center.x > mesh_max_x) {
                mesh_max_x = grid_cell->center.x;
                mesh_max_dx = grid_cell->discretization.x;
            }
            else if(grid_cell->center.x < mesh_min_x) {
                mesh_min_x = grid_cell->center.x;
                mesh_min_dx = grid_cell->discretization.x;
            }

            if(grid_cell->center.y > mesh_max_y) {
                mesh_max_y = grid_cell->center.y;
                mesh_max_dy = grid_cell->discretization.y;
            }
            else if(grid_cell->center.y < mesh_min_y) {
                mesh_min_y = grid_cell->center.y;
                mesh_min_dy = grid_cell->discretization.y;
            }

            if(grid_cell->center.z > mesh_max_z) {
                mesh_max_z = grid_cell->center.z;
                mesh_max_dz = grid_cell->discretization.z;
            }
            else if(grid_cell->center.z < mesh_min_z) {
                mesh_min_z = grid_cell->center.z;
                mesh_min_dz = grid_cell->discretization.z;
            }

        }
    }

    result.x = mesh_max_x;
    if(mesh_max_x != mesh_min_x)
        result.x = (mesh_max_x-mesh_min_x)/2.0f;


    result.y = mesh_max_y;
    if(mesh_max_y != mesh_min_y)
        result.y = (mesh_max_y-mesh_min_y)/2.0f;

    result.z = mesh_max_z;
    if(mesh_max_z != mesh_min_z)
        result.z = (mesh_max_z-mesh_min_z)/2.0f;

    mesh_info->center_calculated = true;

    mesh_info->max_size.x = mesh_max_x + mesh_max_dx/2.0f;
    mesh_info->max_size.y = mesh_max_y + mesh_max_dy/2.0f;
    mesh_info->max_size.z = mesh_max_z + mesh_max_dz/2.0f;

    mesh_info->min_size.x = mesh_min_x - mesh_min_dx/2.0f;
    mesh_info->min_size.y = mesh_min_y - mesh_min_dy/2.0f;
    mesh_info->min_size.z = mesh_min_z - mesh_min_dz/2.0f;

    return result;

}

static Vector3 find_mesh_center_vtk(struct mesh_info *mesh_info) {

    struct vtk_unstructured_grid *grid_to_draw = gui_config.grid_info.vtk_grid;

    uint32_t n_active = grid_to_draw->num_cells;

    float mesh_max_x, mesh_max_y, mesh_max_z;
    float mesh_min_x, mesh_min_y, mesh_min_z;

    mesh_max_x = FLT_MIN;
    mesh_max_y = FLT_MIN;
    mesh_max_z = FLT_MIN;

    mesh_min_x = FLT_MAX;
    mesh_min_y = FLT_MAX;
    mesh_min_z = FLT_MAX;

    int64_t *cells = grid_to_draw->cells;
    point3d_array points = grid_to_draw->points;

    Vector3 result = (Vector3){0.0f, 0.0f, 0.0f};
    float dx, dy, dz;
    float mesh_center_x, mesh_center_y, mesh_center_z;
    uint32_t num_points = grid_to_draw->points_per_cell;

    float max_dx, max_dy, max_dz;
    float min_dx, min_dy, min_dz;

    max_dx = FLT_MIN;
    max_dy = FLT_MIN;
    max_dz = FLT_MIN;

    min_dx = FLT_MAX;
    min_dy = FLT_MAX;
    min_dz = FLT_MAX;

    for (uint32_t i = 0; i < n_active*num_points; i+=num_points) {

        dx = fabs((points[cells[i]].x - points[cells[i+1]].x));
        dy = fabs((points[cells[i]].y - points[cells[i+3]].y));
        dz = fabs((points[cells[i]].z - points[cells[i+4]].z));

        mesh_center_x = points[cells[i]].x + dx/2.0f;
        mesh_center_y = points[cells[i]].y + dy/2.0f;
        mesh_center_z = points[cells[i]].z + dz/2.0f;

        if(mesh_center_x > mesh_max_x) {
            mesh_max_x = mesh_center_x;
            max_dx = dx;
        }
        else if(mesh_center_x < mesh_min_x) {
            mesh_min_x = mesh_center_x;
            min_dx = dx;
        }

        if(mesh_center_y > mesh_max_y) {
            mesh_max_y = mesh_center_y;
            max_dy = dy;
        }
        else if(mesh_center_y < mesh_min_y) {
            mesh_min_y = mesh_center_y;
            min_dy = dy;
        }

        if(mesh_center_z > mesh_max_z) {
            mesh_max_z = mesh_center_z;
            max_dz = dz;
        }
        else if(mesh_center_z < mesh_min_z) {
            mesh_min_z = mesh_center_z;
            min_dz = dz;
        }

    }

    result.x = mesh_max_x;
    if(mesh_max_x != mesh_min_x)
        result.x = (mesh_max_x-mesh_min_x)/2.0f;


    result.y = mesh_max_y;
    if(mesh_max_y != mesh_min_y)
        result.y = (mesh_max_y-mesh_min_y)/2.0f;

    result.z = mesh_max_z;
    if(mesh_max_z != mesh_min_z)
        result.z = (mesh_max_z-mesh_min_z)/2.0f;

    mesh_info->center_calculated = true;

    mesh_info->max_size.x = mesh_max_x + max_dx/2.0f;
    mesh_info->max_size.y = mesh_max_y + max_dy/2.0f;
    mesh_info->max_size.z = mesh_max_z + max_dz/2.0f;

    mesh_info->min_size.x = mesh_min_x - min_dx/2.0f;
    mesh_info->min_size.y = mesh_min_y - min_dy/2.0f;
    mesh_info->min_size.z = mesh_min_z - min_dz/2.0f;

    return result;

}

static void draw_voxel(struct voxel *voxel, struct gui_state *gui_state) {

    Vector3 p_draw = voxel->position_draw;
    Vector3 p_mesh = voxel->position_mesh;

	Vector3 cube_size = voxel->size;
	real_cpu v = voxel->v;

    bool collision = false;
    bool collision_mouse_over = false;

    real_cpu max_v = gui_config.max_v;
    real_cpu min_v = gui_config.min_v;

    Color color;

    action_potential_array aps = (struct action_potential*) hmget(gui_state->ap_graph_config->selected_aps, p_draw);

    struct action_potential ap1;
    ap1.t = gui_config.time;
    ap1.v = voxel->v;
    size_t aps_len = arrlen(aps);

    if(aps != NULL) {
        if(ap1.t > aps[aps_len-1].t ) {
            arrput(aps, ap1);
            hmput(gui_state->ap_graph_config->selected_aps, p_draw, aps);
        }
    }

	if(gui_state->double_clicked) {
		collision = CheckCollisionRayBox( gui_state->ray,
				(BoundingBox){(Vector3){p_draw.x - cube_size.x / 2, p_draw.y - cube_size.y / 2, p_draw.z - cube_size.z / 2},
				(Vector3){p_draw.x + cube_size.x / 2, p_draw.y + cube_size.y / 2, p_draw.z + cube_size.z / 2}});
	}

	collision_mouse_over = CheckCollisionRayBox(
				gui_state->ray_mouse_over, (BoundingBox){(Vector3){p_draw.x - cube_size.x / 2, p_draw.y - cube_size.y / 2, p_draw.z - cube_size.z / 2},
				(Vector3){p_draw.x + cube_size.x / 2, p_draw.y + cube_size.y / 2, p_draw.z + cube_size.z / 2}});

	if(collision_mouse_over) {
		gui_state->current_mouse_over_volume.position_draw = (Vector3){p_mesh.x, p_mesh.y, p_mesh.z};
		gui_state->current_mouse_over_volume.matrix_position = voxel->matrix_position;
	}


    if(collision && !gui_state->one_selected) {
        gui_state->current_selected_volume = (Vector3){p_mesh.x, p_mesh.y, p_mesh.z};
        gui_state->one_selected = true;
        gui_state->ap_graph_config->draw_selected_ap_text = true;
    }

    color = get_color((v - min_v)/(max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);

    if(gui_state->draw_grid_only) {
        DrawCubeWiresV(p_draw, cube_size, color);
    }
    else {

        DrawCubeV(p_draw, cube_size, color);

        if(gui_state->draw_grid_lines) {
            DrawCubeWiresV(p_draw, cube_size, BLACK);
        }

        if(gui_state->current_selected_volume.x == p_mesh.x && gui_state->current_selected_volume.y == p_mesh.y && gui_state->current_selected_volume.z == p_mesh.z) {

			DrawCubeWiresV(p_draw, cube_size, GREEN);

            if(aps == NULL) {
                arrsetcap(aps, 50);
                struct action_potential ap;
                ap.t = gui_config.time;
                ap.v = v;

                arrput(aps, ap);
                hmput(gui_state->ap_graph_config->selected_aps, p_draw, aps);
                gui_state->selected_time = GetTime();
            }
        }
    }
}

static void draw_vtk_unstructured_grid(Vector3 mesh_offset, real_cpu scale, struct gui_state *gui_state) {

    struct vtk_unstructured_grid *grid_to_draw = gui_config.grid_info.vtk_grid;

    if(!grid_to_draw) return;

    int64_t *cells = grid_to_draw->cells;
    point3d_array points = grid_to_draw->points;

    uint32_t n_active = grid_to_draw->num_cells;

    int num_points = grid_to_draw->points_per_cell;
    int j = num_points;
	
	struct voxel voxel;

    for (uint32_t i = 0; i < n_active*num_points; i+=num_points) {

        float mesh_center_x, mesh_center_y, mesh_center_z;
        float dx, dy, dz;

        dx = fabs((points[cells[i]].x - points[cells[i+1]].x));
        dy = fabs((points[cells[i]].y - points[cells[i+3]].y));
        dz = fabs((points[cells[i]].z - points[cells[i+4]].z));

        mesh_center_x = points[cells[i]].x + dx/2.0f;
        mesh_center_y = points[cells[i]].y + dy/2.0f;
        mesh_center_z = points[cells[i]].z + dz/2.0f;

        voxel.v = grid_to_draw->values[j-num_points];
        j += 1;

        voxel.position_draw.x = (float)((mesh_center_x - mesh_offset.x)/scale);
        voxel.position_draw.y = (float)((mesh_center_y - mesh_offset.y)/scale);
        voxel.position_draw.z = (float)((mesh_center_z - mesh_offset.z)/scale);

        voxel.size.x = (float)(dx/scale);
        voxel.size.y = (float)(dy/scale);
        voxel.size.z = (float)(dz/scale);

		voxel.position_mesh = (Vector3){mesh_center_x, mesh_center_y, mesh_center_z};

        draw_voxel(&voxel, gui_state);
    }
    gui_state->one_selected = false;
}

static void draw_alg_mesh(Vector3 mesh_offset, real_cpu scale, struct gui_state *gui_state) {

    struct grid *grid_to_draw = gui_config.grid_info.alg_grid;

    if(!grid_to_draw) return;

	struct voxel voxel;

    if (grid_to_draw) {

        uint32_t n_active = grid_to_draw->num_active_cells;
        struct cell_node **ac = grid_to_draw->active_cells;

        if (ac) {
            for (uint32_t i = 0; i < n_active; i++) {

                struct cell_node *grid_cell;

                grid_cell = ac[i];

                //if(!cell_is_visible(grid_cell)) {
                //    continue;
               // }

                voxel.position_draw.x = (float)((grid_cell->center.x - mesh_offset.x)/scale);
                voxel.position_draw.y = (float)((grid_cell->center.y - mesh_offset.y)/scale);
                voxel.position_draw.z = (float)((grid_cell->center.z - mesh_offset.z)/scale);

                voxel.size.x = (float)(grid_cell->discretization.x/scale);
                voxel.size.y = (float)(grid_cell->discretization.y/scale);
                voxel.size.z = (float)(grid_cell->discretization.z/scale);
				
				voxel.position_mesh = (Vector3){grid_cell->center.x, grid_cell->center.y, grid_cell->center.z};
				voxel.v = ac[i]->v;
				voxel.matrix_position = ac[i]->grid_position;
                draw_voxel(&voxel, gui_state);

            }
        }
    }
    gui_state->one_selected = false;
}

static inline double clamp(double x, double min, double max) {
    if (x < min)
        x = min;
    else if (x > max)
        x = max;
    return x;
}

static void draw_ap_graph(struct gui_state *gui_state, Font font) {

    static const Color colors[] = {DARKGRAY, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BROWN, DARKBROWN, BLACK, MAGENTA};
    int num_colors = SIZEOF(colors);

    float spacing_big = (float)gui_state->font_size_big/(float)font.baseSize;
    float spacing_small = (float)gui_state->font_size_small/(float)font.baseSize;

    DrawRectangleRec(gui_state->ap_graph_config->graph, WHITE);

    gui_state->ap_graph_config->drag_graph_button_position = (Rectangle){(float)(gui_state->ap_graph_config->graph.x + gui_state->ap_graph_config->graph.width)-7.5f, (float)(gui_state->ap_graph_config->graph.y)-7.5f, 15.0f, 15.0f};
    gui_state->ap_graph_config->move_graph_button_position = (Rectangle){(float)(gui_state->ap_graph_config->graph.x), (float)(gui_state->ap_graph_config->graph.y) - 7.5f, 15.0f, 15.0f};

    GuiButton(gui_state->ap_graph_config->drag_graph_button_position, " ");
    GuiButton(gui_state->ap_graph_config->move_graph_button_position, "+");

    char tmp[1024];
    sprintf(tmp, "%.2lf", gui_config.final_time);

    Vector2 text_width = MeasureTextEx(font, tmp, (float)gui_state->font_size_small, spacing_small);
    gui_state->ap_graph_config->min_x = (float)gui_state->ap_graph_config->graph.x + 2.5f*text_width.x;
    gui_state->ap_graph_config->max_x = (float)gui_state->ap_graph_config->graph.x + (float)gui_state->ap_graph_config->graph.width - text_width.x;

    gui_state->ap_graph_config->min_y = (float) gui_state->ap_graph_config->graph.y + (float)gui_state->ap_graph_config->graph.height - 50.0f;
    gui_state->ap_graph_config->max_y = (float) gui_state->ap_graph_config->graph.y + 20.0f; //This is actually the smallest allowed y

    int n = hmlen(gui_state->ap_graph_config->selected_aps);

    if(gui_state->ap_graph_config->draw_selected_ap_text) {
        char *ap_text = "%d AP(s) selected ( cell at %f, %f, %f )";
        double time_elapsed = GetTime() - gui_state->selected_time;
        unsigned char alpha = (unsigned char) clamp(255 - time_elapsed*25, 0, 255);
        Color c = colors[(n-1) % num_colors];
        c.a = alpha;

        sprintf(tmp, ap_text, n, gui_state->current_selected_volume.x, gui_state->current_selected_volume.y, gui_state->current_selected_volume.z);

        text_width = MeasureTextEx(font, ap_text, (float)gui_state->font_size_big, spacing_big);

        Vector2 text_position = (Vector2){gui_state->ap_graph_config->graph.x + gui_state->ap_graph_config->graph.width/2 - text_width.x/2,
                                          gui_state->ap_graph_config->graph.y + 5};

        DrawTextEx(font, tmp, text_position, gui_state->font_size_big, spacing_big, c);

        if(alpha == 0) {
            gui_state->ap_graph_config->draw_selected_ap_text = false;
            gui_state->selected_time = 0.0;
        }
    }

    //Draw y label
    {
        char *label;
        label = "Vm (mV)";
        text_width = MeasureTextEx(font, label, gui_state->font_size_big, spacing_big);

        Vector2 text_position = (Vector2) {
                gui_state->ap_graph_config->graph.x + 15,
                (float) gui_state->ap_graph_config->min_y - (float)((gui_state->ap_graph_config->min_y - gui_state->ap_graph_config->max_y)/2.0f) + (text_width.x/2.0)};

        DrawTextEx2(font, label, text_position, gui_state->font_size_big, spacing_big, BLACK);
    }

    Vector2 p1, p2;

    uint num_ticks;
    real_cpu tick_ofsset = 10;
    num_ticks = (uint) ((gui_config.max_v - gui_config.min_v)/tick_ofsset);

    if(num_ticks < MIN_VERTICAL_TICKS) {
        num_ticks = MIN_VERTICAL_TICKS;
        tick_ofsset = (gui_config.max_v - gui_config.min_v)/num_ticks;
    }
    else if(num_ticks > MAX_VERTICAL_TICKS) {
        num_ticks = MAX_VERTICAL_TICKS;
        tick_ofsset = (gui_config.max_v - gui_config.min_v)/num_ticks;
    }

    real_cpu v = gui_config.min_v;
    sprintf(tmp, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(font, tmp, (float)gui_state->font_size_small, spacing_small);

    //Draw vertical ticks (Vm)
    for(uint t = 0; t <= num_ticks; t++ ) {

        p1.x = (float)gui_state->ap_graph_config->graph.x + 5.0f;
        p1.y = normalize(gui_config.min_v, gui_config.max_v, gui_state->ap_graph_config->min_y, gui_state->ap_graph_config->max_y, v);

        sprintf(tmp, "%.2lf",  normalize(gui_state->ap_graph_config->min_y, gui_state->ap_graph_config->max_y, gui_config.min_v, gui_config.max_v, p1.y));
        text_width = MeasureTextEx(font, tmp, (float)gui_state->font_size_small, spacing_small);

        DrawTextEx(font, tmp, (Vector2){p1.x + (max_w.x - text_width.x/2)  + 20, p1.y - text_width.y/2.0f}, (float)gui_state->font_size_small, spacing_small, RED);

        p1.x = gui_state->ap_graph_config->min_x - 5.0f;
        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        //Grid line
        DrawLineV(p1, (Vector2){gui_state->ap_graph_config->max_x, p1.y}, LIGHTGRAY);

        //Tick line
        DrawLineV(p1, p2, RED);

        v += tick_ofsset;
    }

    tick_ofsset = 10;
    num_ticks = (int) (gui_config.final_time/tick_ofsset);

    if(num_ticks < MIN_HORIZONTAL_TICKS) {
        num_ticks = MIN_HORIZONTAL_TICKS;
        tick_ofsset = gui_config.final_time/num_ticks;
    }
    else if(num_ticks > MAX_HORIZONTAL_TICKS) {
        num_ticks = MAX_HORIZONTAL_TICKS;
        tick_ofsset = gui_config.final_time/num_ticks;
    }

    real_cpu time = 0.0;

    //Draw horizontal ticks (t)
    for(uint t = 0; t <= num_ticks; t++) {

        p1.x = normalize(0.0f, gui_config.final_time, gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, time);
        p1.y = gui_state->ap_graph_config->min_y - 5;

        p2.x = p1.x;
        p2.y = gui_state->ap_graph_config->min_y + 5;

        if(!(t%2)) {
            sprintf(tmp, "%.2lf", normalize(gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, 0.0f, gui_config.final_time, p1.x));
            text_width = MeasureTextEx(font, tmp, (float) gui_state->font_size_small, spacing_small);
            DrawTextEx(font, tmp, (Vector2){p1.x - text_width.x/2.0f, p1.y + 10}, (float) gui_state->font_size_small, spacing_small, RED);
        }

        DrawLineV(p1, p2, RED);

        DrawLineV(p1, (Vector2){p1.x, gui_state->ap_graph_config->max_y}, LIGHTGRAY);
        time+= tick_ofsset;
    }

    //Draw x label
    {
        char *time_text;

        if (gui_config.dt == 0) {
            time_text = "Time (steps)";
        } else {
            time_text = "Time (ms)";
        }
        text_width = MeasureTextEx(font, time_text, gui_state->font_size_big, spacing_big);

        Vector2 text_position = (Vector2) {
                gui_state->ap_graph_config->graph.x + (float) gui_state->ap_graph_config->graph.width / 2.0f -
                text_width.x / 2.0f,
                (float) gui_state->ap_graph_config->min_y + 25.0f};

        DrawTextEx(font, time_text, text_position, gui_state->font_size_big, spacing_big, BLACK);
    }

    //Draw vertical line
    {
        p1.x = gui_state->ap_graph_config->min_x;
        p1.y = gui_state->ap_graph_config->min_y + 5;

        p2.x = p1.x;
        p2.y = gui_state->ap_graph_config->max_y;
        DrawLineV(p1, p2, RED);
    }

    //Draw horizontal line
    {
        p1.x = gui_state->ap_graph_config->min_x;
        p1.y = gui_state->ap_graph_config->min_y;

        p2.x = gui_state->ap_graph_config->max_x;
        p2.y = p1.y;
        DrawLineV(p1, p2, RED);
    }

    struct action_potential *aps;

    //Draw the function
    for (int j = 0; j < n; j++) {

        aps = (struct action_potential*) gui_state->ap_graph_config->selected_aps[j].value;
        int c = arrlen(aps);
        int step = 1;

        if(c > 0) {
            Color line_color = colors[j % num_colors];
            for (int i = 0; i < c; i+=step) {

                if(aps[i].t <= gui_config.time) {
                    if (i + step < c) {

                        p1.x = normalize(0.0f, gui_config.final_time, gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, aps[i].t);
                        p1.y = normalize(gui_config.min_v, gui_config.max_v, gui_state->ap_graph_config->min_y, gui_state->ap_graph_config->max_y, aps[i].v);

                        p2.x = normalize(0.0f, gui_config.final_time, gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, aps[i + step].t);
                        p2.y = normalize(gui_config.min_v, gui_config.max_v, gui_state->ap_graph_config->min_y, gui_state->ap_graph_config->max_y, aps[i + step].v);

                        //TODO: create an option for this???
                        if (aps[i + step].v > gui_config.max_v) gui_config.max_v = aps[i + step].v;
                        if (aps[i + step].v < gui_config.min_v) gui_config.min_v = aps[i + step].v;

                        DrawLineV(p1, p2, line_color);
                    }
                }

            }
        }
    }

    //Draw AP coordinates over mouse cursor
    if(!gui_state->ap_graph_config->drag_ap_graph && gui_state->ap_graph_config->selected_ap_point.x != FLT_MAX && gui_state->ap_graph_config->selected_ap_point.y != FLT_MAX) {
        char *tmp_point = "%lf, %lf";
        sprintf(tmp, tmp_point, gui_state->ap_graph_config->selected_ap_point.x, gui_state->ap_graph_config->selected_ap_point.y);
        text_width = MeasureTextEx(font, tmp, (float)gui_state->font_size_small, spacing_big);
        DrawTextEx(font, tmp, (Vector2){gui_state->mouse_pos.x-text_width.x/2, gui_state->mouse_pos.y-text_width.y}, (float)gui_state->font_size_small, 1, BLACK);
    }

    if(gui_state->ap_graph_config->selected_point_for_apd1.x != FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd1.y != FLT_MAX) {
        DrawCircleV(gui_state->ap_graph_config->selected_point_for_apd1, 4, RED);
    }

    if(gui_state->ap_graph_config->selected_point_for_apd2.x != FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd2.y != FLT_MAX) {
        DrawCircleV(gui_state->ap_graph_config->selected_point_for_apd2, 4, RED);
        DrawLineV(gui_state->ap_graph_config->selected_point_for_apd1, gui_state->ap_graph_config->selected_point_for_apd2, RED);

        float t1 = normalize(gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, 0.0f, gui_config.final_time, gui_state->ap_graph_config->selected_point_for_apd1.x);
        float t2 = normalize(gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, 0.0f, gui_config.final_time, gui_state->ap_graph_config->selected_point_for_apd2.x);

        char *tmp_point = "dt = %lf";
        sprintf(tmp, tmp_point, fabsf(t2-t1));
        text_width = MeasureTextEx(font, tmp, (float)gui_state->font_size_small, spacing_big);
        DrawTextEx(font, tmp, (Vector2){gui_state->ap_graph_config->selected_point_for_apd1.x+text_width.x/2.0, gui_state->ap_graph_config->selected_point_for_apd1.y-text_width.y}, (float)gui_state->font_size_small, 1, BLACK);
    }
}

static inline void move_rect(Vector2 new_pos, Rectangle *rect) {

    float new_x = new_pos.x;
    float new_y = new_pos.y;

    if(new_x > 1 && new_x + rect->width < GetScreenWidth()) {
        rect->x = new_x;
    }
    if(new_y > 1 && new_y + rect->height < GetScreenHeight()) {
        rect->y = new_y;
    }
}

static inline void drag_box(Vector2 mouse_pos, Rectangle *box) {
    float new_x = mouse_pos.x - box->width / 2;
    float new_y = mouse_pos.y - box->height / 2;
    move_rect((Vector2){new_x, new_y}, box);
}

static inline void drag_scale(Vector2 new_pos, Rectangle *box) {
    float new_x = new_pos.x - box->width / 2;
    float new_y = new_pos.y + box->height / 2;
    move_rect((Vector2){new_x, new_y}, box);
}


static void check_window_bounds(Rectangle *box) {
    if (box->x + box->width > current_window_width)
        move_rect((Vector2){current_window_width - box->width - 10, box->y}, box);

    if (box->y + box->height > current_window_height)
        move_rect((Vector2){box->x, current_window_height - box->height - 10}, box);
}

static void draw_scale(struct gui_state *gui_state, bool int_scale) {

    static const int scale_width = 20;
    check_window_bounds(&(gui_state->scale_bounds));

    static bool calc_bounds = true;

    float spacing_small = (float)gui_state->font_size_small/(float)gui_state->font.baseSize;
    float spacing_big = (float)gui_state->font_size_big/(float)gui_state->font.baseSize;

    int num_ticks;
    real_cpu tick_ofsset = 12;

    if(!int_scale) {
        num_ticks = (int) ((gui_config.max_v - gui_config.min_v) / tick_ofsset);

        if (num_ticks < 5) {
            num_ticks = 5;
            tick_ofsset = (gui_config.max_v - gui_config.min_v) / num_ticks;
        }
    }
    else {
        num_ticks = 0;
        for(int i = gui_config.min_v; i < gui_config.max_v; i++) {
            num_ticks++;
        }
        tick_ofsset = (gui_config.max_v - gui_config.min_v) / num_ticks;
    }

    char tmp[256];

    real_cpu v = gui_config.min_v;
    sprintf(tmp, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(gui_state->font, tmp, gui_state->font_size_small, spacing_small);

    Vector2 p1, p2, width;

    float scale_rec_height  = 30.0f;
    Color color;

    real_cpu  min_v = gui_config.min_v;
    real_cpu  max_v = gui_config.max_v;

    if(calc_bounds) {
        gui_state->scale_bounds.height += max_w.y;
        for(int t = 0; t <= num_ticks; t++ ) {
            gui_state->scale_bounds.height += scale_rec_height;
        }
        calc_bounds = false;
    }

    if(!int_scale) {
        width = MeasureTextEx(gui_state->font, "Vm", gui_state->font_size_big, spacing_big);
        float diff = scale_width - width.x;

        p1.x = gui_state->scale_bounds.x + (diff/2.0);
        p1.y = (float)gui_state->scale_bounds.y - (float)gui_state->scale_bounds.height + 20;
        DrawTextEx(gui_state->font, "Vm", p1, gui_state->font_size_big, spacing_big, BLACK);
    }

    float initial_y = gui_state->scale_bounds.y;

    for(int t = 0; t <= num_ticks; t++ ) {
        p1.x = gui_state->scale_bounds.x - 55.0f;
        p1.y = initial_y + (float)scale_rec_height/2.0f;

        sprintf(tmp, "%.2lf", v);
        width = MeasureTextEx(gui_state->font, tmp, gui_state->font_size_small, spacing_small);

        DrawTextEx(gui_state->font, tmp, (Vector2){p1.x + (max_w.x - width.x), p1.y - width.y/2.0f}, gui_state->font_size_small, spacing_small, BLACK);

        p1.x = p1.x + max_w.x + 2.5f;
        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        DrawLineV(p1, p2, BLACK);
        color = get_color((v - min_v)/(max_v - min_v), gui_state->scale_alpha, gui_state->current_scale);

        DrawRectangle((int)gui_state->scale_bounds.x,(int)initial_y, scale_width, (int)scale_rec_height, color);
        initial_y -= scale_rec_height;
        v += tick_ofsset;
    }

}

static void draw_box(Rectangle *box, int text_offset, const char **lines, int num_lines, int font_size_for_line, Font font) {

    check_window_bounds(box);

    int text_x = (int)box->x + 20;
    int text_y = (int)box->y + 10;

    DrawRectangleRec(*box, WHITE);

    DrawRectangleLinesEx(*box, 1, BLACK);

    for (int i = 0; i < num_lines; i++) {
        DrawTextEx(font, lines[i], (Vector2){text_x, text_y}, font_size_for_line, 1, BLACK);
        text_y += text_offset;
    }
}

static inline void configure_end_info_box_strings (char ***info_string) {

    char tmp[128];

    int index = 0;

    sprintf(tmp, "Resolution Time: %ld s", gui_config.solver_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "ODE Total Time: %ld s", gui_config.ode_total_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "CG Total Time: %ld s", gui_config.cg_total_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Mat time: %ld s", gui_config.total_mat_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Refine time: %ld s", gui_config.total_ref_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Unrefine time: %ld s", gui_config.total_deref_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Write time: %ld s", gui_config.total_write_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Initial configuration time: %ld s", gui_config.total_config_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "CG Total Iterations: %ld", gui_config.total_cg_it);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Final Time: %lf ms", gui_config.time);
    (*(info_string))[index]  = strdup(tmp);

}

static inline bool configure_mesh_info_box_strings (char ***info_string, int draw_type, struct mesh_info *mesh_info) {

    if(!gui_config.grid_info.alg_grid && !gui_config.grid_info.vtk_grid) return false;

    char tmp[128];

    int index = 0;

    uint32_t n_active = 0;

    if(draw_type == DRAW_SIMULATION) {
        n_active = gui_config.grid_info.alg_grid->num_active_cells;
    }
    else {
        n_active = gui_config.grid_info.vtk_grid->num_cells;
    }

    (*(info_string))[index++] = strdup("Mesh information:");

    sprintf(tmp, " - Num. of Volumes: %u", n_active);
    (*(info_string))[index++] = strdup(tmp);

    sprintf(tmp, " - Max X: %f", mesh_info->max_size.x);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Max Y: %f", mesh_info->max_size.y);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Max Z: %f", mesh_info->max_size.z);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Min X: %f", mesh_info->min_size.x);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Min Y: %f", mesh_info->min_size.y);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Min Z: %f", mesh_info->min_size.z);
    (*(info_string))[index++]  = strdup(tmp);

    if(draw_type == DRAW_SIMULATION) {
        if (gui_config.paused) {
            sprintf(tmp, "Simulation paused: %lf of %lf ms", gui_config.time, gui_config.final_time);
        } else if (gui_config.simulating) {
            sprintf(tmp, "Simulation running: %lf of %lf ms", gui_config.time, gui_config.final_time);

        } else {
            sprintf(tmp, "Simulation finished: %lf of %lf ms", gui_config.time, gui_config.final_time);
        }
    }
    else {
        if (gui_config.paused) {
            if(gui_config.dt == 0) {
                sprintf(tmp, "Visualization paused: %d of %d steps", (int)gui_config.time, (int)gui_config.final_time);
            }
            else {
                sprintf(tmp, "Visualization paused: %lf of %lf ms", gui_config.time, gui_config.final_time);
            }
        } else if (gui_config.simulating) {
            if(gui_config.dt == 0) {
                sprintf(tmp, "Visualization running: %d of %d steps", (int) gui_config.time, (int) gui_config.final_time);
            }
            else {
                sprintf(tmp, "Visualization running: %lf of %lf ms", gui_config.time, gui_config.final_time);
            }

        } else {
            if(gui_config.dt == 0) {
                sprintf(tmp, "Visualization finished: %d of %d steps", (int)gui_config.time, (int)gui_config.final_time);
            }
            else {
                sprintf(tmp, "Visualization finished: %lf of %lf ms", gui_config.time, gui_config.final_time);
            }
        }
    }

    (*(info_string))[index] = strdup(tmp);

    return true;
}

static bool draw_selection_box(struct gui_state *gui_state) {

    const float text_box_width = 60;
    const float text_box_height = 25;
    const float text_box_y_dist = 40;
    const float label_box_y_dist = 30;
    const float x_off = 10;

    static char center_x_text[128] = { 0 };
    static char center_y_text[128] = { 0 };
    static char center_z_text[128] = { 0 };

    float pos_x = gui_state->sub_window_pos.x;
    float pos_y = gui_state->sub_window_pos.y;

    float box_pos = pos_x + x_off;

    bool window_closed = GuiWindowBox((Rectangle){ pos_x, pos_y, gui_state->box_width , gui_state->box_height}, "Enter the center of the cell");

    DrawTextEx(gui_state->font, "Center X", (Vector2){box_pos + 5, pos_y + label_box_y_dist}, gui_state->font_size_small, 1, BLACK);
    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_width, text_box_height}, center_x_text, SIZEOF(center_x_text) - 1, true);

    box_pos = pos_x + text_box_width + 2*x_off;
    DrawTextEx(gui_state->font, "Center Y", (Vector2){box_pos + 5, pos_y + label_box_y_dist}, gui_state->font_size_small, 1, BLACK);
    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_width, text_box_height}, center_y_text, SIZEOF(center_y_text) - 1, true);

    box_pos = pos_x +  2*text_box_width + 3*x_off;
    DrawTextEx(gui_state->font, "Center Z", (Vector2){box_pos + 5, pos_y + label_box_y_dist}, gui_state->font_size_small, 1, BLACK);
    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_width, text_box_height}, center_z_text, SIZEOF(center_z_text) - 1, true);

    bool btn_ok_clicked = GuiButton((Rectangle){pos_x + text_box_width + 2*x_off, pos_y + 70, text_box_width, text_box_height}, "OK");

    if(btn_ok_clicked) {
        gui_state->current_selected_volume.x = strtof(center_x_text, NULL);
        gui_state->current_selected_volume.y = strtof(center_y_text, NULL);
        gui_state->current_selected_volume.z = strtof(center_z_text, NULL);
    }

    return window_closed || btn_ok_clicked;
}

static void reset(struct mesh_info *mesh_info, struct gui_state *gui_state, bool full_reset) {

    gui_state->voxel_alpha = 255;

    for(long  i = 0; i < hmlen(gui_state->ap_graph_config->selected_aps); i++) {
        arrsetlen(gui_state->ap_graph_config->selected_aps[i].value, 0);
    }

    if(full_reset) {
        for(long i = 0; i < hmlen(gui_state->ap_graph_config->selected_aps); i++) {
            arrfree(gui_state->ap_graph_config->selected_aps[i].value);
        }
        hmfree(gui_state->ap_graph_config->selected_aps);
        gui_state->ap_graph_config->selected_aps = NULL;
        hmdefault(gui_state->ap_graph_config->selected_aps, NULL);

        set_camera_params(&(gui_state->camera));
    }

    gui_state->ap_graph_config->draw_selected_ap_text = false;

    if(gui_config.paused) {
        omp_unset_lock(&gui_config.sleep_lock);
        gui_config.paused = false;
    }

    gui_config.restart = true;
    gui_config.grid_info.alg_grid = NULL;
    gui_config.grid_info.vtk_grid = NULL;

    gui_state->ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    omp_unset_lock(&gui_config.sleep_lock);
    gui_state->current_selected_volume = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->current_mouse_over_volume.position_draw = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd1 = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd2 = (Vector2){FLT_MAX, FLT_MAX};
}

static void handle_keyboard_input(struct mesh_info *mesh_info, struct gui_state *gui_state) {

    if(gui_config.paused) {

        if(IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown((KEY_LEFT_CONTROL))) {
            //SAVE FILE AS VTK
            if(IsKeyPressed(KEY_S)) {
                char const * filter[1] = {"*.vtu"};

                const char *save_path = tinyfd_saveFileDialog ("Save VTK file", gui_config.input, 1, filter, "vtu files" );

                if(save_path) {
                    save_vtk_unstructured_grid_as_vtu_compressed(gui_config.grid_info.vtk_grid, save_path, 6);
                    log_to_stdout_and_file("Saved vtk file as %s\n", save_path);
                }

                return;
            }
        }

        if (IsKeyPressed(KEY_RIGHT) || IsKeyDown(KEY_UP)) {
            gui_config.advance_or_return = 1;
            omp_unset_lock(&gui_config.sleep_lock);
            nanosleep((const struct timespec[]){{0, 50000000L}}, NULL);
            return;
        }

        if(gui_config.draw_type == DRAW_FILE) {
            //Return one step only works on file visualization...
            if (IsKeyPressed(KEY_LEFT) || IsKeyDown(KEY_DOWN)) {
                gui_config.advance_or_return = -1;
                nanosleep((const struct timespec[]){{0, 50000000L}}, NULL);
                omp_unset_lock(&gui_config.sleep_lock);
                return;
            }
        }

    }

    if(IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown((KEY_LEFT_CONTROL))) {
        if(IsKeyPressed(KEY_F)) {
            gui_state->show_selection_box = true;
            gui_state->sub_window_pos.x = GetScreenWidth() / 2 - gui_state->box_width;
            gui_state->sub_window_pos.y = GetScreenHeight() / 2 - gui_state->box_height;
            return;
        }
    }

    if (IsKeyPressed(KEY_Q)) {
        gui_state->show_scale = !gui_state->show_scale;
        return;
    }

    if (IsKeyDown(KEY_A)) {
        if(gui_state->voxel_alpha - 1 >= 0) {
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

    if (IsKeyPressed(KEY_G))  {
        gui_state->draw_grid_only = !gui_state->draw_grid_only;
        return;
    }

    if (IsKeyPressed(KEY_PERIOD))  {
        gui_state->current_scale = (gui_state->current_scale + 1) % NUM_SCALES;
        return;
    }

    if (IsKeyPressed(KEY_COMMA))  {
        if(gui_state->current_scale - 1 >= 0) {
            gui_state->current_scale = (gui_state->current_scale - 1);
        }
        else {
            gui_state->current_scale = NUM_SCALES-1;
        }
        return;
    }

    if (IsKeyPressed(KEY_L)) {
        gui_state->draw_grid_lines = !gui_state->draw_grid_lines;
        return;
    }

    if (IsKeyPressed(KEY_SPACE)) {
        gui_config.paused = !gui_config.paused;
        return;
    }

    if (IsKeyPressed(KEY_R)) {
        bool full_reset = false;

        if(IsKeyDown(KEY_LEFT_ALT)) {
            full_reset = true;
        }

        reset(mesh_info, gui_state, full_reset);
        return;
    }

    if (IsKeyPressed(KEY_X)) {
        gui_state->show_ap = !gui_state->show_ap;
        return;
    }

    if (IsKeyPressed(KEY_C)) {
        gui_state->show_scale = gui_state->c_pressed;
        gui_state->show_ap = gui_state->c_pressed;
        gui_state->show_help_box = gui_state->c_pressed;
        gui_state->show_end_info_box = gui_state->c_pressed;
        gui_state->show_mesh_info_box = gui_state->c_pressed;
        gui_state->c_pressed = !gui_state->c_pressed;
        gui_state->show_coordinates = !gui_state->show_coordinates;
        return;
    }

    if (IsKeyPressed(KEY_O)) {

        gui_config.paused = true;

        char *buf = get_current_directory();

        char const *tmp = tinyfd_selectFolderDialog("Select a directory", buf);
        if(tmp) {
            gui_config.input = strdup(tmp);
        }
        else {
            gui_config.input = NULL;
        }

        free(buf);

        if(tmp) {
            reset(mesh_info, gui_state, true);
        }

        return;
    }

    if (IsKeyPressed(KEY_F)) {

        gui_config.paused = true;

        char *buf = get_current_directory();

        char const * filter[4] = {"*.pvd", "*.acm", "*.vtk", "*.vtu"};

        char const *tmp = tinyfd_openFileDialog(
            "Select a simulation file",
            buf,
            4,
            filter,
            "simulation result (pvd, vtk, vtu or acm)",
            0);

        if(tmp) {
            gui_config.input = strdup(tmp);
        }
        else {
            gui_config.input = NULL;
        }

        free(buf);

        if(tmp) {
            reset(mesh_info, gui_state, true);
        }
        return;
    }

}

static void handle_input(struct mesh_info *mesh_info, struct gui_state *gui_state) {

    if(gui_state->handle_keyboard_input) {
        handle_keyboard_input(mesh_info, gui_state);
    }

    gui_state->mouse_pos = GetMousePosition();
	gui_state->ray_mouse_over = GetMouseRay(GetMousePosition(), gui_state->camera);

    if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {

    	gui_state->ray = GetMouseRay(GetMousePosition(), gui_state->camera);

        if(!gui_state->show_selection_box) {

            if (gui_state->mouse_timer == -1) {

                gui_state->double_clicked = false;
                gui_state->mouse_timer = GetTime();

            } else {

                double delay = GetTime() - gui_state->mouse_timer;

                if (delay < DOUBLE_CLICK_DELAY) {
                    gui_state->double_clicked = true;
                    gui_state->mouse_timer = -1;
                } else {
                    gui_state->mouse_timer = -1;
                    gui_state->double_clicked = false;
                }
            }
        }

        if (CheckCollisionRayBox(gui_state->ray,
                                 (BoundingBox){(Vector3){ gui_state->coordinates_cube.x - gui_state->coordinates_cube_size.x/2.0,
                                                          gui_state->coordinates_cube.y - gui_state->coordinates_cube_size.y/2.0,
                                                          gui_state->coordinates_cube.z - gui_state->coordinates_cube_size.z/2.0},
                                               (Vector3){ gui_state->coordinates_cube.x + gui_state->coordinates_cube_size.x/2.0,
                                                          gui_state->coordinates_cube.y + gui_state->coordinates_cube_size.y/2.0,
                                                          gui_state->coordinates_cube.z + gui_state->coordinates_cube_size.z/2.0}})) {
            gui_state->move_coordinates = true;
        }

        else if (CheckCollisionPointRec(gui_state->mouse_pos, (Rectangle){gui_state->sub_window_pos.x, gui_state->sub_window_pos.y,
                                                                     gui_state->box_width - 18, WINDOW_STATUSBAR_HEIGHT })) {
            gui_state->move_sub_window = true;
        }
        else if (CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->drag_graph_button_position)) {
            gui_state->ap_graph_config->drag_ap_graph = true;
        }
        else if (CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->move_graph_button_position)) {
            gui_state->ap_graph_config->move_ap_graph = true;
        }
        else if (CheckCollisionPointRec(gui_state->mouse_pos, gui_state->help_box)) {
            gui_state->move_help_box = true;
        }
        else if (CheckCollisionPointRec(gui_state->mouse_pos, gui_state->mesh_info_box)) {
            gui_state->move_info_box = true;
        }
        else if (CheckCollisionPointRec(gui_state->mouse_pos, gui_state->end_info_box)) {
            gui_state->move_end_info_box = true;
        }
        else if (CheckCollisionPointRec(gui_state->mouse_pos,
                (Rectangle){gui_state->scale_bounds.x, gui_state->scale_bounds.y - gui_state->scale_bounds.height,
                            gui_state->scale_bounds.width, gui_state->scale_bounds.height})) {
            gui_state->move_scale = true;
        }

    }

    else if(IsMouseButtonPressed(MOUSE_RIGHT_BUTTON)) {
        if(hmlen(gui_state->ap_graph_config->selected_aps) && gui_state->show_ap) {
            if (CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->graph)) {
                if(gui_state->ap_graph_config->selected_point_for_apd1.x == FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd1.y == FLT_MAX) {
                    gui_state->ap_graph_config->selected_point_for_apd1.x = gui_state->mouse_pos.x;
                    gui_state->ap_graph_config->selected_point_for_apd1.y = gui_state->mouse_pos.y;
                }
                else {
                    if(gui_state->ap_graph_config->selected_point_for_apd2.x == FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd2.y == FLT_MAX) {
                        gui_state->ap_graph_config->selected_point_for_apd2.x = gui_state->mouse_pos.x;
                        gui_state->ap_graph_config->selected_point_for_apd2.y = gui_state->ap_graph_config->selected_point_for_apd1.y;
                    }
                    else {
                        gui_state->ap_graph_config->selected_point_for_apd1.x = gui_state->mouse_pos.x;
                        gui_state->ap_graph_config->selected_point_for_apd1.y = gui_state->mouse_pos.y;

                        gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
                        gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;
                    }
                }
            }

        }
    }

    if (gui_state->move_sub_window) {
        gui_state->sub_window_pos.x = (gui_state->mouse_pos.x) - (gui_state->box_width - 18) / 2;
        gui_state->sub_window_pos.y = (gui_state->mouse_pos.y) - WINDOW_STATUSBAR_HEIGHT / 2;
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) gui_state->move_sub_window = false;
    }

    else if (gui_state->ap_graph_config->drag_ap_graph) {

        float new_heigth  = gui_state->ap_graph_config->graph.height + (gui_state->ap_graph_config->graph.y -  gui_state->mouse_pos.y);

        if(new_heigth > 100){
            gui_state->ap_graph_config->graph.height = new_heigth;
            gui_state->ap_graph_config->graph.y = gui_state->mouse_pos.y;
        }

        float new_width = gui_state->mouse_pos.x - gui_state->ap_graph_config->graph.x;

        if(new_width > 200) {
            gui_state->ap_graph_config->graph.width = new_width;
        }

        gui_state->ap_graph_config->selected_point_for_apd1.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd1.y = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;

        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) gui_state->ap_graph_config->drag_ap_graph = false;
    }

    else if (gui_state->ap_graph_config->move_ap_graph) {

        if(gui_state->mouse_pos.y > 10 && gui_state->mouse_pos.x + gui_state->ap_graph_config->graph.width < (float)GetScreenWidth()) {
            gui_state->ap_graph_config->graph.x = gui_state->mouse_pos.x;
        }

        if(gui_state->mouse_pos.y > 10 && gui_state->mouse_pos.y + gui_state->ap_graph_config->graph.height < (float)GetScreenHeight()) {
            gui_state->ap_graph_config->graph.y = gui_state->mouse_pos.y;
        }

        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) gui_state->ap_graph_config->move_ap_graph = false;

        gui_state->ap_graph_config->selected_point_for_apd1.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd1.y = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;
    }

    else if (gui_state->move_help_box) {
        drag_box(gui_state->mouse_pos, &gui_state->help_box);
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) gui_state->move_help_box = false;
    }
    else if(gui_state->move_info_box) {
        drag_box(gui_state->mouse_pos, &gui_state->mesh_info_box);
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) gui_state->move_info_box = false;
    }
    else if(gui_state->move_scale) {
        drag_scale(gui_state->mouse_pos, &gui_state->scale_bounds);
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) gui_state->move_scale = false;
    }
    else if(gui_state->move_end_info_box) {
        drag_box(gui_state->mouse_pos, &gui_state->end_info_box);
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) gui_state->move_end_info_box = false;

    }

    if(hmlen(gui_state->ap_graph_config->selected_aps) && gui_state->show_ap) {
        float t = normalize(gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, 0.0f, gui_config.final_time, gui_state->mouse_pos.x);
        float v = normalize(gui_state->ap_graph_config->min_y, gui_state->ap_graph_config->max_y, gui_config.min_v, gui_config.max_v, gui_state->mouse_pos.y);
        if (CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->graph)) {
            gui_state->ap_graph_config->selected_ap_point.x = t;
            gui_state->ap_graph_config->selected_ap_point.y = v;
        }
        else {
            gui_state->ap_graph_config->selected_ap_point.x = FLT_MAX;
            gui_state->ap_graph_config->selected_ap_point.y = FLT_MAX;
        }

    }

}

static int configure_info_boxes_sizes(struct gui_state *gui_state) {
    Vector2 txt_w_h;
    int text_offset;
    int box_w = 0;

    txt_w_h = MeasureTextV(WIDER_TEXT, (int)gui_state->font_size_small);
    text_offset = (int)(1.5 * txt_w_h.y);
    box_w = (int)(txt_w_h.x + 50);

    gui_state->help_box.width = (float)box_w;
    gui_state->help_box.height = (float)(text_offset * info_box_lines) + 10.0f;

    gui_state->mesh_info_box.width = (float)box_w;
    gui_state->mesh_info_box.height = (float)(text_offset * mesh_info_box_lines) + 10;

    gui_state->mesh_info_box.x = (float)(current_window_width - box_w - 10);
    gui_state->mesh_info_box.y = 10.0f;

    gui_state->end_info_box.width = (float)box_w;
    gui_state->end_info_box.height = (float)(text_offset * end_info_box_lines) + 10;
    gui_state->end_info_box.x = gui_state->help_box.x;
    gui_state->end_info_box.y = gui_state->help_box.y + gui_state->help_box.height + 10;

    return text_offset;
}

void draw_coordinates(struct gui_state *gui_state) {

    const float line_size = 1.0f;
    const float arrow_offset = 0.1f;
    static bool first_draw = true;

    if(first_draw) {
        gui_state->coordinates_cube = (Vector3){-(line_size / 2.0) + 0.5, -7 + 0.5, -1.5};
        first_draw = false;
    }

    Vector3 start_pos = (Vector3){gui_state->coordinates_cube.x-0.5, gui_state->coordinates_cube.y-0.5, gui_state->coordinates_cube.z-0.5};
    Vector3 end_pos = (Vector3){start_pos.x + line_size, start_pos.y, gui_state->coordinates_cube.z-0.5};

    DrawLine3D(start_pos, end_pos,  RED);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y + arrow_offset, end_pos.z}, end_pos, RED);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, RED);

    gui_state->coordinates_label_x_position =  GetWorldToScreen(end_pos, gui_state->camera);

    end_pos = (Vector3){start_pos.x, start_pos.y + line_size, end_pos.z};

    DrawLine3D(start_pos, end_pos,  GREEN);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, GREEN);
    DrawLine3D((Vector3){end_pos.x + arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, GREEN);

    gui_state->coordinates_label_y_position =  GetWorldToScreen((Vector3){end_pos.x, end_pos.y + 0.2, end_pos.z}, gui_state->camera);

    end_pos = (Vector3){start_pos.x, start_pos.y, start_pos.z + line_size};

    DrawLine3D(start_pos, end_pos,  DARKBLUE);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y, end_pos.z - arrow_offset}, end_pos, DARKBLUE);
    DrawLine3D((Vector3){end_pos.x + arrow_offset, end_pos.y, end_pos.z - arrow_offset}, end_pos, DARKBLUE);

    gui_state->coordinates_label_z_position =  GetWorldToScreen(end_pos, gui_state->camera);
    //DrawCubeWiresV(gui_state->coordinates_cube, (Vector3){1.2,1.2,1.2}, WHITE);

}


void init_and_open_gui_window() {

    omp_set_lock(&gui_config.sleep_lock);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT);
    char *window_title = NULL;

    int draw_type = gui_config.draw_type;

    if(draw_type == DRAW_SIMULATION) {
        window_title = (char *) malloc(strlen(gui_config.config_name) + strlen("Simulation visualization - ") + 2);
        sprintf(window_title, "Simulation visualization - %s", gui_config.config_name);
    }
    else {
        window_title = strdup("Opening mesh...");
    }

    int current_monitor = 0;


    InitWindow(current_window_width, current_window_height, window_title);

    current_monitor = GetCurrentMonitor();
    current_window_width = GetMonitorWidth(current_monitor);
    current_window_height = GetMonitorHeight(current_monitor);

    SetWindowSize(current_window_width, current_window_height);

    free(window_title);

    SetTargetFPS(60);

    Image icon = LoadImage("res/icon.png");

    if(icon.data)
        SetWindowIcon(icon);

    UnloadImage(icon);

    float scale = 1.0f;

    Vector3 mesh_offset = (Vector3){ 0, 0, 0 };

    char **end_info_box_strings = NULL;
    char **mesh_info_box_strings = NULL;

    const char *info_box_strings[] = {
            "Default controls:",
            " - Mouse Wheel to Zoom in-out",
            " - Mouse Wheel Pressed to Pan",
            " - Alt + Mouse Wheel Pressed to Rotate",
            " - Alt + Ctrl + Mouse Wheel Pressed for Smooth Zoom",
            " - Ctrl + F to search a cell based on it's center",
            " - G to only draw the grid lines",
            " - L to enable or disable the grid lines",
            " - R to restart simulation",
            " - X to show/hide AP visualization",
            " - Q to show/hide scale",
            " - C to show/hide everything except grid",
            " - F to open a simulation file",
            " - O to open a simulation directory",
            " - . or , to change color scales",
            " - Right arrow to advance one dt when paused",
            " - Hold up arrow to advance time when paused",
            " - Double click on a volume to show the AP",
            " - Space to start or pause simulation"
    };

    info_box_lines = SIZEOF(info_box_strings);

    end_info_box_strings = (char **) malloc(sizeof(char *) * end_info_box_lines);
    mesh_info_box_strings = (char **) malloc(sizeof(char *) * mesh_info_box_lines);

    Vector2 error_message_witdh;

    struct mesh_info *mesh_info = new_mesh_info();
    struct gui_state *gui_state = new_gui_state_with_font_sizes(14, 18);

    bool end_info_box_strings_configured = false;

    int text_offset = configure_info_boxes_sizes(gui_state);

    while (!WindowShouldClose()) {

        UpdateCamera(&(gui_state->camera));
		


        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

        if (IsWindowResized()) {
            current_window_height = GetScreenHeight();
            current_window_width = GetScreenWidth();
        }

        gui_state->handle_keyboard_input = !gui_state->show_selection_box;


        handle_input(mesh_info, gui_state);
        if(gui_config.grid_info.loaded) {

            omp_set_lock(&gui_config.draw_lock);

            if(draw_type == DRAW_FILE) {
                if(gui_config.grid_info.file_name) {
                    window_title = (char *) malloc(strlen(gui_config.grid_info.file_name) + strlen("Visualizing file - ") + 2);
                    sprintf(window_title, "Visualizing file - %s", gui_config.grid_info.file_name);
                    SetWindowTitle(window_title);
                    free(window_title);
                }
            }

            ClearBackground(GRAY);

            BeginMode3D(gui_state->camera);

            if(!mesh_info->center_calculated) {

                if(draw_type == DRAW_SIMULATION) {
                    mesh_offset = find_mesh_center(mesh_info);
					 //scale = fmaxf(gui_config.grid_info.alg_grid->mesh_side_length.x,
                     //             fmaxf(gui_config.grid_info.alg_grid->mesh_side_length.y,
                     //                   gui_config.grid_info.alg_grid->mesh_side_length.z)) / 5.0f;
                }
                else {
                    mesh_offset = find_mesh_center_vtk(mesh_info);
                }
				
				scale = fmaxf(mesh_offset.x, fmaxf(mesh_offset.y, mesh_offset.z)) / 5.0f;
					

            }

            if(draw_type == DRAW_SIMULATION) {
                draw_alg_mesh(mesh_offset, scale, gui_state);
            }
            else if(draw_type == DRAW_FILE) {
                draw_vtk_unstructured_grid(mesh_offset, scale, gui_state);
            }

            if(gui_state->show_coordinates) {
                draw_coordinates(gui_state);
            }

            EndMode3D();

            if(gui_state->show_coordinates) {
                DrawText("x", gui_state->coordinates_label_x_position.x, gui_state->coordinates_label_x_position.y, gui_state->font_size_big, RED);
                DrawText("y", gui_state->coordinates_label_y_position.x, gui_state->coordinates_label_y_position.y, gui_state->font_size_big, GREEN);
                DrawText("z", gui_state->coordinates_label_z_position.x, gui_state->coordinates_label_z_position.y, gui_state->font_size_big, DARKBLUE);
            }

            if(gui_state->show_mesh_info_box) {
                bool configured = configure_mesh_info_box_strings(&mesh_info_box_strings, draw_type, mesh_info);

                if(configured) {
                    draw_box(&gui_state->mesh_info_box, text_offset, (const char **) mesh_info_box_strings,
                             mesh_info_box_lines, (int)gui_state->font_size_small, gui_state->font);

                    for (int i = 0; i < mesh_info_box_lines; i++) {
                        free(mesh_info_box_strings[i]);
                    }
                }
            }
            //We finished drawing everything that depends on the mesh being loaded
            omp_unset_lock(&gui_config.draw_lock);

            if(gui_state->show_scale) {
                draw_scale(gui_state, gui_config.int_scale);
            }

            if(hmlen(gui_state->ap_graph_config->selected_aps) && gui_state->show_ap) {
                draw_ap_graph(gui_state, gui_state->font);
            }

           if(gui_state->show_help_box) {
               draw_box(&gui_state->help_box, text_offset, info_box_strings, info_box_lines, (int)gui_state->font_size_small, gui_state->font);
           }

            if(!gui_config.simulating) {
                if(draw_type == DRAW_SIMULATION) {
                    if (gui_state->show_end_info_box) {
                        if(!end_info_box_strings_configured) {
                            configure_end_info_box_strings(&end_info_box_strings);
                            end_info_box_strings_configured = true;
                        }
                        draw_box(&gui_state->end_info_box, text_offset, (const char **) end_info_box_strings,
                                 end_info_box_lines, (int)gui_state->font_size_small, gui_state->font);
                    }
                }
            }

            if(!gui_config.paused) {
                omp_unset_lock(&gui_config.sleep_lock);
            }

            if(gui_state->show_selection_box) {
                gui_state->show_selection_box = !draw_selection_box(gui_state);
            }

        }
        else {

            ClearBackground(GRAY);
            float spacing = 20/(float)gui_state->font.baseSize;
            static Color c = RED;

            if(!gui_config.error_message) {
                gui_config.error_message = strdup("Loading Mesh...");
                c = WHITE;
            }

            error_message_witdh = MeasureTextEx(gui_state->font, gui_config.error_message, 20, spacing);

            int posx = GetScreenWidth()/2 - (int)error_message_witdh.x/2;
            int posy = GetScreenHeight()/2 - 50;

            int rec_width = (int)(error_message_witdh.x) + 40;

            DrawRectangle(posx, posy, rec_width, 20, c);
            DrawRectangleLines(posx, posy, rec_width, 20, BLACK);
            if(gui_config.error_message) //This should not happen... but it does....
                DrawText(gui_config.error_message, posx + 20, posy, 20, BLACK);

        }

        //Draw FPS
        int fps = GetFPS();
        DrawText(TextFormat("%2i FPS", fps), GetScreenWidth()  - 100, GetScreenHeight()-20, 20, BLACK);
        DrawText(TextFormat("Mouse is on Volume: %lf, %lf, %lf with grid position %i", 
				gui_state->current_mouse_over_volume.position_draw.x, gui_state->current_mouse_over_volume.position_draw.y, 
				gui_state->current_mouse_over_volume.position_draw.z, gui_state->current_mouse_over_volume.matrix_position), 
				GetScreenWidth() - MeasureText("Mouse is on Volume: 10000, 10000, 10000 with grid position 10", 20)*1.4, GetScreenHeight()-50, 20, BLACK);
        EndDrawing();

    }

    gui_config.exit = true;
    free(mesh_info);

    if(end_info_box_strings_configured) {
        for (int i = 0; i < end_info_box_lines; i++) {
            free(end_info_box_strings[i]);
        }
    }

    free(end_info_box_strings);

    hmfree(gui_state->ap_graph_config->selected_aps);

    free(gui_state->ap_graph_config);
    free(gui_state);

    omp_unset_lock(&gui_config.draw_lock);
    omp_unset_lock(&gui_config.sleep_lock);

    CloseWindow();

}
