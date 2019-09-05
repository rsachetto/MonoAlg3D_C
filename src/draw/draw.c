//
// Created by sachetto on 11/11/17.
//

#include <float.h>
#include <time.h>
#include <limits.h>

#include "draw.h"
#include "color_maps.h"
#include "../single_file_libraries/stb_ds.h"
#include "../raylib/src/raylib.h"
#include "../raylib/src/camera.h"
#include "../raylib/src/raymath.h"
#include "../raylib/src/rlgl.h"

#define RAYGUI_IMPLEMENTATION
#define RAYGUI_RICONS_SUPPORT
#define RAYGUI_TEXTBOX_EXTENDED
#include "../raylib/src/raygui.h"

static bool calc_center = false;
static bool one_selected = false;
static bool draw_selected_ap_text = false;
static bool show_ap = true;
static bool show_scale = true;
static bool c_pressed = false;
static bool draw_grid_lines = false;
static bool draw_grid_only = false;

static bool show_info_box = true;
static bool show_end_info_box = true;
static bool show_mesh_info_box = true;
static bool show_selection_box = false;
static bool show_save_box = false;

static Vector3 max_size;
static Vector3 min_size;

#define SIZEOF(A) (sizeof(A)/sizeof(A[0]))
#define WIDER_TEXT "------------------------------------------------------"
#define DOUBLE_CLICK_DELAY 0.5 //seconds

static double mouse_timer = -1;

static double selected_time = 0.0;
static Vector3 current_selected = {FLT_MAX, FLT_MAX, FLT_MAX};

char center_x_text[128] = { 0 };
char center_y_text[128] = { 0 };
char center_z_text[128] = { 0 };
char save_path[PATH_MAX] = {0};

static float center_x;
static float center_y;
static float center_z;

static Vector2 mouse_pos = {0 };
static Vector2 selected_ap_point = { FLT_MAX, FLT_MAX };
static Vector2 selected_point_for_apd1 = { FLT_MAX, FLT_MAX };
static Vector2 selected_point_for_apd2 = { FLT_MAX, FLT_MAX };
static Vector2 window_pos;
static bool drag_window = false;

//AP VARIABLES
static int graph_height;
static int graph_pos_x;
static int graph_pos_y;
static int graph_width;

static float max_x;
static float min_x;
static float min_y;
static float max_y;

static int current_scale = 0;

struct point_voidp_hash_entry *selected_aps;

static Color colors[] = {DARKGRAY, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BROWN, DARKBROWN, BLACK, MAGENTA};
static int num_colors = SIZEOF(colors);

static inline float normalize(float r_min, float r_max, float t_min, float t_max, float m) {
    return ((m - r_min) / (r_max-r_min))*(t_max - t_min) + t_min;
}

static inline Color get_color(real_cpu value)
{

    int idx1;        // |-- Our desired color will be between these two indexes in "color".
    int idx2;        // |
    real_cpu fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

    if(value <= 0)      {  idx1 = idx2 = 0;            }    // accounts for an input <=0
    else if(value >= 1)  {  idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
    else
    {
        value = value * (NUM_COLORS-1);        // Will multiply value by 3.
        idx1  = (int)floor(value);                  // Our desiBLACK color will be after this index.
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
    result.a = 255;

    return result;
}

static inline bool skip_node(struct cell_node *grid_cell) {

    if(!cell_has_neighbour(grid_cell, grid_cell->north) ) {
        return false;
    }
    else if(!cell_has_neighbour(grid_cell, grid_cell->south) ) {
        return false;
    }
    else if(!cell_has_neighbour(grid_cell, grid_cell->west) ) {
        return false;
    }
    else if(!cell_has_neighbour(grid_cell, grid_cell->east) ) {
        return false;
    }
    else if(!cell_has_neighbour(grid_cell, grid_cell->front) ) {
        return false;
    }
    else if(!cell_has_neighbour(grid_cell, grid_cell->back) ) {
        return false;
    }
    else {
        return true;
    }
}

static Vector3 find_mesh_center() {

    struct grid *grid_to_draw = draw_config.grid_info.grid_to_draw;

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
            if(grid_cell->translated_center.x > mesh_max_x) {
                mesh_max_x = grid_cell->translated_center.x;
                mesh_max_dx = grid_cell->discretization.x;
            }
            else if(grid_cell->translated_center.x < mesh_min_x) {
                mesh_min_x = grid_cell->translated_center.x;
                mesh_min_dx = grid_cell->discretization.x;
            }

            if(grid_cell->translated_center.y > mesh_max_y) {
                mesh_max_y = grid_cell->translated_center.y;
                mesh_max_dy = grid_cell->discretization.y;
            }
            else if(grid_cell->translated_center.y < mesh_min_y) {
                mesh_min_y = grid_cell->translated_center.y;
                mesh_min_dy = grid_cell->discretization.y;
            }

            if(grid_cell->translated_center.z > mesh_max_z) {
                mesh_max_z = grid_cell->translated_center.z;
                mesh_max_dz = grid_cell->discretization.z;
            }
            else if(grid_cell->translated_center.z < mesh_min_z) {
                mesh_min_z = grid_cell->translated_center.z;
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

    calc_center = true;

    max_size.x = mesh_max_x + mesh_max_dx/2.0f;
    max_size.y = mesh_max_y + mesh_max_dy/2.0f;
    max_size.z = mesh_max_z + mesh_max_dz/2.0f;

    min_size.x = mesh_min_x - mesh_min_dx/2.0f;
    min_size.y = mesh_min_y - mesh_min_dy/2.0f;
    min_size.z = mesh_min_z - mesh_min_dz/2.0f;

    return result;

}

static Vector3 find_mesh_center_vtk() {

    struct vtk_unstructured_grid *grid_to_draw = draw_config.grid_info.vtk_grid;

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
    uint32_t num_points =grid_to_draw->points_per_cell;

    float max_dx, max_dy, max_dz;
    float min_dx, min_dy, min_dz;

    max_dx = FLT_MIN;
    max_dy = FLT_MIN;
    max_dz = FLT_MIN;

    min_dx = FLT_MAX;
    min_dy = FLT_MAX;
    min_dz = FLT_MAX;

    for (uint32_t i = 0; i < n_active*num_points; i+=num_points) {

        dx = fabsf((points[cells[i]].x - points[cells[i+1]].x));
        dy = fabsf((points[cells[i]].y - points[cells[i+3]].y));
        dz = fabsf((points[cells[i]].z - points[cells[i+4]].z));


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

    calc_center = true;

    max_size.x = mesh_max_x + max_dx/2.0f;
    max_size.y = mesh_max_y + max_dy/2.0f;
    max_size.z = mesh_max_z + max_dz/2.0f;

    min_size.x = mesh_min_x - min_dx/2.0f;
    min_size.y = mesh_min_y - min_dy/2.0f;
    min_size.z = mesh_min_z - min_dz/2.0f;

    return result;

}

static void draw_voxel(Vector3 cube_position_draw, Vector3 cube_position_mesh, Vector3 cube_size, real_cpu v, Ray ray) {

    struct point_3d p;
    p.x = cube_position_draw.x;
    p.y = cube_position_draw.y;
    p.z = cube_position_draw.z;

    bool collision;

    real_cpu max_v = draw_config.max_v;
    real_cpu min_v = draw_config.min_v;

    Color color;

    action_potential_array aps = (struct action_potential*) hmget(selected_aps, p);

    struct action_potential ap1;
    ap1.t = draw_config.time;
    ap1.v = v;
    size_t aps_len = arrlen(aps);

    if(aps != NULL) {
        if(ap1.t > aps[aps_len-1].t ) {
            arrput(aps, ap1);
            hmput(selected_aps, p, aps);
        }
    }

    collision = CheckCollisionRayBox(ray,
                                     (BoundingBox){(Vector3){ cube_position_draw.x - cube_size.x/2, cube_position_draw.y - cube_size.y/2, cube_position_draw.z - cube_size.z/2 },
                                                   (Vector3){ cube_position_draw.x + cube_size.x/2, cube_position_draw.y + cube_size.y/2, cube_position_draw.z + cube_size.z/2 }});

    if(collision && !one_selected) {
        current_selected = (Vector3){cube_position_mesh.x, cube_position_mesh.y, cube_position_mesh.z};
        one_selected = true;
        draw_selected_ap_text = true;
    }

    color = get_color((v - min_v)/(max_v - min_v));

    if(draw_grid_only) {
        DrawCubeWiresV(cube_position_draw, cube_size, color);
    }
    else {

        DrawCubeV(cube_position_draw, cube_size, color);

        if(draw_grid_lines) {
            DrawCubeWiresV(cube_position_draw, cube_size, BLACK);
        }

        if(current_selected.x == cube_position_mesh.x && current_selected.y == cube_position_mesh.y && current_selected.z == cube_position_mesh.z) {
            DrawCubeWiresV(cube_position_draw, cube_size, GREEN);

            if(aps == NULL) {
                arrsetcap(aps, 50);
                struct action_potential ap;
                ap.t = draw_config.time;
                ap.v = v;

                arrput(aps, ap);
                hmput(selected_aps, p, aps);
                selected_time = GetTime();
            }
        }
    }
}

static void draw_vtk_unstructured_grid(Vector3 mesh_offset, real_cpu scale, Ray ray) {

    struct vtk_unstructured_grid *grid_to_draw = draw_config.grid_info.vtk_grid;


    Vector3 cube_position;
    Vector3 cube_size;


    int64_t *cells = grid_to_draw->cells;
    point3d_array points = grid_to_draw->points;

    uint32_t n_active = grid_to_draw->num_cells;

    int num_points = grid_to_draw->points_per_cell;
    int j = num_points;

    for (uint32_t i = 0; i < n_active*num_points; i+=num_points) {

        float mesh_center_x, mesh_center_y, mesh_center_z;
        real_cpu v;
        float dx, dy, dz;

        dx = fabsf((points[cells[i]].x - points[cells[i+1]].x));
        dy = fabsf((points[cells[i]].y - points[cells[i+3]].y));
        dz = fabsf((points[cells[i]].z - points[cells[i+4]].z));

        mesh_center_x = points[cells[i]].x + dx/2.0f;
        mesh_center_y = points[cells[i]].y + dy/2.0f;
        mesh_center_z = points[cells[i]].z + dz/2.0f;

        v = grid_to_draw->values[j-num_points];
        j += 1;

        cube_position.x = (float)((mesh_center_x - mesh_offset.x)/scale);
        cube_position.y = (float)((mesh_center_y - mesh_offset.y)/scale);
        cube_position.z = (float)((mesh_center_z - mesh_offset.z)/scale);

        cube_size.x = (float)(dx/scale);
        cube_size.y = (float)(dy/scale);
        cube_size.z = (float)(dz/scale);

        draw_voxel(cube_position, (Vector3){mesh_center_x, mesh_center_y, mesh_center_z}, cube_size, v, ray);

    }
    one_selected = false;
}

static void draw_alg_mesh(Vector3 mesh_offset, real_cpu scale, Ray ray) {

    struct grid *grid_to_draw = draw_config.grid_info.grid_to_draw;

    Vector3 cube_position;
    Vector3 cube_size;

    if (grid_to_draw) {

        uint32_t n_active = grid_to_draw->num_active_cells;
        struct cell_node **ac = grid_to_draw->active_cells;

        if (ac) {
            for (uint32_t i = 0; i < n_active; i++) {

                struct cell_node *grid_cell;

                grid_cell = ac[i];

                if(skip_node(grid_cell)) {
                    continue;
                }

                cube_position.x = (float)((grid_cell->translated_center.x - mesh_offset.x)/scale);
                cube_position.y = (float)((grid_cell->translated_center.y - mesh_offset.y)/scale);
                cube_position.z = (float)((grid_cell->translated_center.z - mesh_offset.z)/scale);

                cube_size.x = (float)(grid_cell->discretization.x/scale);
                cube_size.y = (float)(grid_cell->discretization.y/scale);
                cube_size.z = (float)(grid_cell->discretization.z/scale);

                draw_voxel(cube_position, (Vector3){grid_cell->translated_center.x, grid_cell->translated_center.y, grid_cell->translated_center.z}, cube_size, ac[i]->v, ray);

            }
        }
    }
    one_selected = false;
}

double clamp(double x, double min, double max) {
    if (x < min)
        x = min;
    else if (x > max)
        x = max;
    return x;
}

void draw_ap(Font font, float font_size_small, float font_size_big) {

    graph_height = GetScreenHeight()/3;

    graph_pos_x = 1;
    graph_pos_y = GetScreenHeight() - graph_height;

    graph_width = GetScreenWidth();

    float spacing_big = (float)font_size_big/(float)font.baseSize;
    float spacing_small = (float)font_size_small/(float)font.baseSize;

    DrawRectangle(graph_pos_x, graph_pos_y, graph_width, graph_height, WHITE);

    char tmp[1024];
    sprintf(tmp, "%.2lf", draw_config.final_time);

    Vector2 width = MeasureTextEx(font, tmp, (float)font_size_small, spacing_small);

    min_x = (float)graph_pos_x + 2.0f*width.x;
    max_x = (float)graph_pos_x + (float)graph_width - width.x;

    min_y = (float) graph_pos_y + (float)graph_height - 30.0f;
    max_y = (float) graph_pos_y + 20.0f; //This is actually the smallest allowed y

    int n = hmlen(selected_aps);

    if(draw_selected_ap_text) {
        char *ap_text = "%d AP(s) selected ( cell at %f, %f, %f )";
        double time_elapsed = GetTime() - selected_time;
        unsigned char alpha = (unsigned char) clamp(255 - time_elapsed*25, 0, 255);
        Color c = colors[(n-1) % num_colors];
        c.a = alpha;

        sprintf(tmp, ap_text, n, current_selected.x, current_selected.y, current_selected.z);

        width = MeasureTextEx(font, ap_text, (float)font_size_big, spacing_big);

        DrawTextEx(font, tmp, (Vector2){(float)graph_pos_x + (float)graph_width/2.0f - width.x/2.0f - min_x, max_y}, font_size_big, 1, c);

        if(alpha == 0) {
            draw_selected_ap_text = false;
            selected_time = 0.0;
        }
    }

    char *time_text;

    if(draw_config.dt == 0) {
        time_text = "Time (steps)";
    }
    else {
        time_text = "Time (ms)";
    }
    width = MeasureTextEx(font, time_text, font_size_big, spacing_big);

    DrawTextEx(font, time_text, (Vector2){min_x + (float)graph_width/2.0f - width.x/2.0f, (float)min_y + 50.0f}, font_size_big, spacing_big, BLACK);

    Vector2 p1, p2;

    uint num_ticks;
    real_cpu tick_ofsset = 10;
    num_ticks = (uint) ((draw_config.max_v - draw_config.min_v)/tick_ofsset);

    if(num_ticks < MIN_VERTICAL_TICKS) {
        num_ticks = MIN_VERTICAL_TICKS;
        tick_ofsset = (draw_config.max_v - draw_config.min_v)/num_ticks;
    }
    else if(num_ticks > MAX_VERTICAL_TICKS) {
        num_ticks = MAX_VERTICAL_TICKS;
        tick_ofsset = (draw_config.max_v - draw_config.min_v)/num_ticks;

    }

    real_cpu v = draw_config.min_v;
    sprintf(tmp, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(font, tmp, (float)font_size_small, spacing_small);

    //Draw vertical ticks (Vm)
    for(uint t = 0; t <= num_ticks; t++ ) {

        p1.x = (float)graph_pos_x + 5.0f;
        p1.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, v);

        sprintf(tmp, "%.2lf",  normalize(min_y, max_y, draw_config.min_v, draw_config.max_v, p1.y));
        width = MeasureTextEx(font, tmp, (float)font_size_small, spacing_small);

        DrawTextEx(font, tmp, (Vector2){p1.x + (max_w.x - width.x), p1.y - width.y/2.0f}, (float)font_size_small, spacing_small, RED);

        p1.x = min_x - 5.0f;
        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        //Grid line
        DrawLineV(p1, (Vector2){max_x, p1.y}, LIGHTGRAY);

        //Tick line
        DrawLineV(p1, p2, RED);

        v += tick_ofsset;
    }

    tick_ofsset = 10;
    num_ticks = (int) (draw_config.final_time/tick_ofsset);

    if(num_ticks < MIN_HORIZONTAL_TICKS) {
        num_ticks = MIN_HORIZONTAL_TICKS;
        tick_ofsset = draw_config.final_time/num_ticks;
    }
    else if(num_ticks > MAX_HORIZONTAL_TICKS) {
        num_ticks = MAX_HORIZONTAL_TICKS;
        tick_ofsset = draw_config.final_time/num_ticks;
    }

    real_cpu time = 0.0;

    //Draw horizontal ticks (t)
    for(uint t = 0; t <= num_ticks; t++) {

        p1.x = normalize(0.0f, draw_config.final_time, min_x, max_x, time);
        p1.y = min_y - 5;

        p2.x = p1.x;
        p2.y = min_y + 5;

        if(!(t%2)) {
            sprintf(tmp, "%.2lf", normalize(min_x, max_x, 0.0f, draw_config.final_time, p1.x));
            width = MeasureTextEx(font, tmp, (float) font_size_small, spacing_small);
            DrawTextEx(font, tmp, (Vector2){p1.x - width.x/2.0f, p1.y + 10}, (float) font_size_small, spacing_small, RED);
        }

        DrawLineV(p1, p2, RED);

        DrawLineV(p1, (Vector2){p1.x, max_y}, LIGHTGRAY);
        time+= tick_ofsset;
    }

    //Draw vertical line
    {
        p1.x = min_x;
        p1.y = min_y + 5;

        p2.x = p1.x;
        p2.y = max_y;
        DrawLineV(p1, p2, RED);
    }

    //Draw horizontal line
    {
        p1.x = min_x;
        p1.y = min_y;

        p2.x = max_x;
        p2.y = p1.y;
        DrawLineV(p1, p2, RED);
    }

    struct action_potential *aps;

    for (int j = 0; j < n; j++) {

        aps = (struct action_potential*) selected_aps[j].value;
        int c = arrlen(aps);

        if(c > 0) {
            Color line_color = colors[j % num_colors];
            for (int i = 0; i < c; i++) {

                if(aps[i].t <= draw_config.time) {
                    if (i + 1 < c) {

                        p1.x = normalize(0.0f, draw_config.final_time, min_x, max_x, aps[i].t);
                        p1.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, aps[i].v);

                        p2.x = normalize(0.0f, draw_config.final_time, min_x, max_x, aps[i + 1].t);
                        p2.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, aps[i + 1].v);

                        //TODO: create an option for this???
                        if (aps[i + 1].v > draw_config.max_v) draw_config.max_v = aps[i + 1].v;
                        if (aps[i + 1].v < draw_config.min_v) draw_config.min_v = aps[i + 1].v;

                        DrawLineV(p1, p2, line_color);
                    }
                }

            }
        }
    }

    if(selected_ap_point.x != FLT_MAX && selected_ap_point.y != FLT_MAX) {
        char *tmp_point = "%lf, %lf";
        sprintf(tmp, tmp_point, selected_ap_point.x, selected_ap_point.y);
        width = MeasureTextEx(font, tmp, (float)font_size_small, spacing_big);
        DrawTextEx(font, tmp, (Vector2){mouse_pos.x-width.x/2, mouse_pos.y-width.y}, (float)font_size_small, 1, BLACK);
    }

    if(selected_point_for_apd1.x != FLT_MAX && selected_point_for_apd1.y != FLT_MAX) {
        DrawCircleV(selected_point_for_apd1, 4, RED);
    }

    if(selected_point_for_apd2.x != FLT_MAX && selected_point_for_apd2.y != FLT_MAX) {
        DrawCircleV(selected_point_for_apd2, 4, RED);
        DrawLineV(selected_point_for_apd1, selected_point_for_apd2, RED);

        float t1 = normalize(min_x, max_x, 0.0f, draw_config.final_time, selected_point_for_apd1.x);
        float t2 = normalize(min_x, max_x, 0.0f, draw_config.final_time, selected_point_for_apd2.x);

        char *tmp_point = "dt = %lf";
        sprintf(tmp, tmp_point, fabsf(t2-t1));
        width = MeasureTextEx(font, tmp, (float)font_size_small, spacing_big);
        DrawTextEx(font, tmp, (Vector2){selected_point_for_apd1.x+width.x/2.0, selected_point_for_apd1.y-width.y}, (float)font_size_small, 1, BLACK);
    }
}

static void draw_scale(Font font, float font_size_small) {

    float initial_y =  (float)GetScreenHeight()/2.0f;

    float spacing_small = (float)font_size_small/(float)font.baseSize;

    int num_ticks;
    real_cpu tick_ofsset = 12;
    num_ticks = (int) ((draw_config.max_v - draw_config.min_v)/tick_ofsset);

    if(num_ticks < 5) {
        num_ticks = 5;
        tick_ofsset = (draw_config.max_v - draw_config.min_v)/num_ticks;
    }

    char tmp[256];

    real_cpu v = draw_config.min_v;
    sprintf(tmp, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    Vector2 p1, p2, width;
    float scale_pos_x = (float)GetScreenWidth()-30.0f;

    float scale_rec_size = 30.0f;
    Color color;

    real_cpu  min_v = draw_config.min_v;
    real_cpu  max_v = draw_config.max_v;

    for(int t = 0; t <= num_ticks; t++ ) {

        p1.x = scale_pos_x - 55.0f;
        p1.y = initial_y + scale_rec_size/2.0f;

        sprintf(tmp, "%.2lf", v);
        width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

        DrawTextEx(font, tmp, (Vector2){p1.x + (max_w.x - width.x), p1.y - width.y/2.0f}, font_size_small, spacing_small, BLACK);

        p1.x = p1.x + max_w.x + 2.5f;
        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        DrawLineV(p1, p2, BLACK);
        color = get_color((v - min_v)/(max_v - min_v));


        DrawRectangle((int)scale_pos_x,(int)initial_y, 20, 30, color);
        initial_y -= scale_rec_size;

        v += tick_ofsset;
    }
}

static void draw_box(int pos_x, int pos_y, int box_w, int box_h, int text_offset, const char **lines, int num_lines, int font_size_for_line) {

    int text_x = pos_x + 20;

    int text_y = pos_y + 10;

    DrawRectangle(pos_x, pos_y, box_w, box_h, WHITE);

    DrawRectangleLines(pos_x, pos_y, box_w, box_h, BLACK);

    for (int i = 0; i < num_lines; i++) {
        DrawText(lines[i], text_x, text_y, font_size_for_line, BLACK);
        text_y += text_offset;
    }

}

static inline void configure_end_info_box_strings (char ***info_string) {

    char tmp[128];

    int index = 0;

    sprintf(tmp, "Resolution Time: %ld s", draw_config.solver_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "ODE Total Time: %ld s", draw_config.ode_total_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "CG Total Time: %ld s", draw_config.cg_total_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Mat time: %ld s", draw_config.total_mat_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Refine time: %ld s", draw_config.total_ref_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Derefine time: %ld s", draw_config.total_deref_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Write time: %ld s", draw_config.total_write_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Initial configuration time: %ld s", draw_config.total_config_time/1000/1000);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "CG Total Iterations: %ld", draw_config.total_cg_it);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Final Time: %lf ms", draw_config.time);
    (*(info_string))[index]  = strdup(tmp);

}

static inline void configure_mesh_info_box_strings (char ***info_string, int draw_type) {

    char tmp[128];

    int index = 0;

    uint32_t n_active = 0;

    if(draw_type == DRAW_SIMULATION) {
        n_active = draw_config.grid_info.grid_to_draw->num_active_cells;
    }
    else {
        n_active = draw_config.grid_info.vtk_grid->num_cells;
    }

    (*(info_string))[index++] = strdup("Mesh information:");

    sprintf(tmp, " - Num. of Volumes: %u", n_active);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Max X: %f", max_size.x);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Max Y: %f", max_size.y);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Max Z: %f", max_size.z);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Min X: %f", min_size.x);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Min Y: %f", min_size.y);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Min Z: %f", min_size.z);
    (*(info_string))[index++]  = strdup(tmp);


    if(draw_type == DRAW_SIMULATION) {
        if (draw_config.paused) {
            sprintf(tmp, "Simulation paused: %lf of %lf ms", draw_config.time, draw_config.final_time);
        } else if (draw_config.simulating) {
            sprintf(tmp, "Simulation running: %lf of %lf ms", draw_config.time, draw_config.final_time);

        } else {
            sprintf(tmp, "Simulation finished: %lf of %lf ms", draw_config.time, draw_config.final_time);
        }
    }
    else {
        if (draw_config.paused) {
            if(draw_config.dt == 0) {
                sprintf(tmp, "Visualization paused: %d of %d steps", (int)draw_config.time, (int)draw_config.final_time);
            }
            else {
                sprintf(tmp, "Visualization paused: %lf of %lf ms", draw_config.time, draw_config.final_time);
            }
        } else if (draw_config.simulating) {
            if(draw_config.dt == 0) {
                sprintf(tmp, "Visualization running: %d of %d steps", (int) draw_config.time, (int) draw_config.final_time);
            }
            else {
                sprintf(tmp, "Visualization running: %lf of %lf ms", draw_config.time, draw_config.final_time);
            }

        } else {
            if(draw_config.dt == 0) {
                sprintf(tmp, "Visualization finished: %d of %d steps", (int)draw_config.time, (int)draw_config.final_time);
            }
            else {
                sprintf(tmp, "Visualization finished: %lf of %lf ms", draw_config.time, draw_config.final_time);
            }
        }
    }

    (*(info_string))[index] = strdup(tmp);


}

const int box_width = 220;
const int box_height = 100;

static bool draw_selection_box(Font font, float font_size) {

    const int text_box_width = 60;
    const int text_box_height = 25;
    const int text_box_y_dist = 40;
    const int label_box_y_dist = 30;

    const int x_off = 10;

    int pos_x = (int) window_pos.x;
    int pos_y = (int) window_pos.y;

    int box_pos = pos_x + x_off;

    bool clicked = GuiWindowBox((Rectangle){ pos_x, pos_y, box_width , box_height}, "Enter the center of the cell");

    DrawTextEx(font, "Center X", (Vector2){box_pos + 5, pos_y + label_box_y_dist}, font_size, 1, BLACK);
    GuiTextBox((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_width, text_box_height}, center_x_text, SIZEOF(center_x_text) - 1, true);
    center_x = atof(center_x_text);

    box_pos = pos_x + text_box_width + 2*x_off;
    DrawTextEx(font, "Center Y", (Vector2){box_pos + 5, pos_y + label_box_y_dist}, font_size, 1, BLACK);
    GuiTextBox((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_width, text_box_height}, center_y_text, SIZEOF(center_y_text) - 1, true);
    center_y = atof(center_y_text);

    box_pos = pos_x +  2*text_box_width + 3*x_off;
    DrawTextEx(font, "Center Z", (Vector2){box_pos + 5, pos_y + label_box_y_dist}, font_size, 1, BLACK);
    GuiTextBox((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_width, text_box_height}, center_z_text, SIZEOF(center_z_text) - 1, true);
    center_z = atof(center_z_text);

    bool btn_clicked = GuiButton((Rectangle){pos_x + text_box_width + 2*x_off, pos_y + 70, text_box_width, text_box_height}, "OK");

    if(btn_clicked) {

        current_selected.x = center_x;
        current_selected.y = center_y;
        current_selected.z = center_z;
    }

    return clicked | btn_clicked;
}

static bool draw_save_box(Font font, float font_size) {

    const int text_box_width = 90;
    const int text_box_height = 25;
    const int text_box_y_dist = 40;

    const int x_off = 10;

    int pos_x = (int) window_pos.x;
    int pos_y = (int) window_pos.y;

    int box_pos = pos_x + x_off;

    bool clicked = GuiWindowBox((Rectangle){ pos_x, pos_y, box_width , box_height}, "Enter the filename");

    GuiTextBox((Rectangle){box_pos, pos_y + text_box_y_dist, box_width - 2*x_off, text_box_height}, save_path, SIZEOF(save_path) - 1, true);

    bool btn_clicked = GuiButton((Rectangle){pos_x +  x_off, pos_y + 70, text_box_width, text_box_height}, "OK");
    bool btn2_clicked = GuiButton((Rectangle){pos_x + box_width - text_box_width - x_off, pos_y + 70, text_box_width, text_box_height}, "CANCEL");

    if(btn_clicked) {
        if( strlen(save_path) > 0 ) {
            save_vtk_unstructured_grid_as_vtu_compressed(draw_config.grid_info.vtk_grid, save_path, 6);
        }
    }

    return btn2_clicked | clicked | btn_clicked;
}

static void handle_input(bool *mesh_loaded, Ray *ray, Camera3D *camera) {

    {
        mouse_pos = GetMousePosition();

        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON))
        {
            if (CheckCollisionPointRec(mouse_pos, (Rectangle){window_pos.x, window_pos.y, box_width - 18 , WINDOW_STATUSBAR_HEIGHT }))
            {
                drag_window = true;
            }
        }

        if (drag_window)
        {
            window_pos.x = (mouse_pos.x) - (box_width - 18) / 2;
            window_pos.y = (mouse_pos.y) - WINDOW_STATUSBAR_HEIGHT / 2;

            if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) drag_window = false;

        }

        if(hmlen(selected_aps) && show_ap) {
            float t = normalize(min_x, max_x, 0.0f, draw_config.final_time, mouse_pos.x);
            float v = normalize(min_y, max_y, draw_config.min_v, draw_config.max_v, mouse_pos.y);
            if (CheckCollisionPointRec(mouse_pos, (Rectangle) {graph_pos_x, graph_pos_y, graph_width, graph_height})) {
                selected_ap_point.x = t;
                selected_ap_point.y = v;
            }
            else {
                selected_ap_point.x = FLT_MAX;
                selected_ap_point.y = FLT_MAX;
            }
        }

        if(IsMouseButtonPressed(MOUSE_RIGHT_BUTTON)) {
            if(hmlen(selected_aps) && show_ap) {
                if (CheckCollisionPointRec(mouse_pos, (Rectangle) {graph_pos_x, graph_pos_y, graph_width, graph_height})) {
                    if(selected_point_for_apd1.x == FLT_MAX && selected_point_for_apd1.y == FLT_MAX) {
                        selected_point_for_apd1.x = mouse_pos.x;
                        selected_point_for_apd1.y = mouse_pos.y;
                    }
                    else {
                        if(selected_point_for_apd2.x == FLT_MAX && selected_point_for_apd2.y == FLT_MAX) {
                            selected_point_for_apd2.x = mouse_pos.x;
                            selected_point_for_apd2.y = selected_point_for_apd1.y;
                        }
                        else {
                            selected_point_for_apd1.x = mouse_pos.x;
                            selected_point_for_apd1.y = mouse_pos.y;

                            selected_point_for_apd2.x = FLT_MAX;
                            selected_point_for_apd2.y = FLT_MAX;
                        }
                    }
                }

            }
        }
    }


    if(IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown((KEY_LEFT_CONTROL))) {
        if(IsKeyPressed(KEY_F)) {
            show_selection_box = true;
            window_pos.x = GetScreenWidth() / 2 - box_width;
            window_pos.y = GetScreenHeight() / 2 - box_height;
        }
    }

    if (IsKeyPressed('Q')) {
        show_scale = !show_scale;
        return;
    }

    if(draw_config.paused) {

        if(IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown((KEY_LEFT_CONTROL))) {
            if(IsKeyPressed(KEY_S)) {
                show_save_box = true;
                window_pos.x = GetScreenWidth() / 2 - box_width;
                window_pos.y = GetScreenHeight() / 2 - box_height;
            }
        }

        if (IsKeyPressed(KEY_RIGHT) || IsKeyDown(KEY_UP)) {
            draw_config.advance_or_return = 1;
            omp_unset_lock(&draw_config.sleep_lock);
            nanosleep((const struct timespec[]){{0, 50000000L}}, NULL);
            return;
        }

        if(draw_config.draw_type == DRAW_FILE) {
            //Return one step only works on file visualization...
            if (IsKeyPressed(KEY_LEFT) || IsKeyDown(KEY_DOWN)) {
                draw_config.advance_or_return = -1;
                nanosleep((const struct timespec[]){{0, 50000000L}}, NULL);
                omp_unset_lock(&draw_config.sleep_lock);
                return;
            }
        }

    }

    if (IsKeyPressed('G'))  {
        draw_grid_only = !draw_grid_only;
        return;
    }

    if (IsKeyPressed('.'))  {
        current_scale = (current_scale + 1) % NUM_SCALES;
        return;
    }
    if (IsKeyPressed(','))  {
        if(current_scale - 1 >= 0) {
            current_scale = (current_scale - 1);
        }
        else {
            current_scale = NUM_SCALES-1;
        }
        return;
    }

    if (IsKeyPressed('L')) {
        draw_grid_lines = !draw_grid_lines;
        return;
    }

    if (IsKeyPressed(KEY_SPACE)) {
        draw_config.paused = !draw_config.paused;
        return;
    }

    if (IsKeyPressed('R')) {
        for(long  i = 0; i < hmlen(selected_aps); i++) {
            arrfree(selected_aps[i].value);
        }

        hmfree(selected_aps);
        selected_aps = NULL;
        hmdefault(selected_aps, NULL);

        if(draw_config.paused) {
            omp_unset_lock(&draw_config.sleep_lock);
            draw_config.paused = false;
        }

        draw_config.restart = true;
        draw_config.grid_info.grid_to_draw = NULL;
        draw_config.grid_info.vtk_grid = NULL;
        *mesh_loaded = false;

        ray->position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
        ray->direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
        calc_center = false;
        omp_unset_lock(&draw_config.sleep_lock);
        current_selected = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
        selected_point_for_apd1 = (Vector2){FLT_MAX, FLT_MAX};
        selected_point_for_apd2 = (Vector2){FLT_MAX, FLT_MAX};
        return;
    }

    if (IsKeyPressed('A')) {
        show_ap = !show_ap;
        return;
    }

    if (IsKeyPressed('C')) {
        show_scale = c_pressed;
        show_ap = c_pressed;
        show_info_box = c_pressed;
        show_end_info_box = c_pressed;
        show_mesh_info_box = c_pressed;
        c_pressed = !c_pressed;
        return;
    }

    if(!show_selection_box) {
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            if (mouse_timer == -1) {
                mouse_timer = GetTime();
            } else {
                double delay = GetTime() - mouse_timer;

                if (delay < DOUBLE_CLICK_DELAY) {
                    *ray = GetMouseRay(GetMousePosition(), *camera);
                    mouse_timer = -1;
                } else {
                    mouse_timer = -1;
                }

            }
        }
    }
}

void init_and_open_visualization_window() {

    omp_set_lock(&draw_config.sleep_lock);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    char *window_title = NULL;

    int draw_type = draw_config.draw_type;

    if(draw_type == DRAW_SIMULATION) {
        window_title = (char *) malloc(strlen(draw_config.config_name) + strlen("Simulation visualization - ") + 2);
        sprintf(window_title, "Simulation visualization - %s", draw_config.config_name);
    }
    else {
        window_title = strdup("Opening mesh...");
    }

    InitWindow(0, 0, window_title);

    free(window_title);

    int current_monitor = 0;

    Font font = GetFontDefault();

    selected_aps = NULL;
    hmdefault(selected_aps, NULL);

    Camera3D camera;

    camera.position = (Vector3){11.082402f, 6.763101f, 8.921088f};  // Camera position
    camera.target = (Vector3){ 1.081274f, -1.581945f, 0.196326f};
    camera.up       = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy     = 45.0f;                                // Camera field-of-view Y
    camera.type     = CAMERA_PERSPECTIVE;                  // Camera mode type

    SetCameraMode(camera, CAMERA_FREE); // Set a free camera mode

    SetTargetFPS(120);

    float scale = 1.0f;

    Ray ray;
    ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};

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
            " - A to show/hide AP visualization",
            " - Q to show/hide scale",
            " - C to show/hide everything except grid",
            " - . or , to change color scales",
            " - Right arrow to advance one dt when paused",
            " - Hold up arrow to advance time when paused",
            " - Double click on a volume to show the AP",
            " - Space to start or pause simulation"
    };

    int info_box_lines = SIZEOF(info_box_strings);
    int end_info_box_lines = 10;
    int mesh_info_box_lines = 9;

    Vector2 txt_w_h;

    int text_offset;

    int box_w, box_h, info_box_h = 0;

    int pos_x;
    float font_size_small;
    float font_size_big;

    end_info_box_strings = (char **) malloc(sizeof(char *) * end_info_box_lines);
    mesh_info_box_strings = (char **) malloc(sizeof(char *) * mesh_info_box_lines);

    Vector2 error_message_witdh;

    while (!WindowShouldClose()) {
        //Configure font size according to monitor resolution
        current_monitor = GetCurrentMonitor();

        font_size_small  = (10.0f*((float)GetMonitorHeight(current_monitor)/1080.0f));
        if(font_size_small < 10) font_size_small = 10;

        font_size_big    = (16.0f*((float)GetMonitorHeight(current_monitor)/1080.0f));
        if(font_size_big < 16) font_size_big = 16;

        txt_w_h = MeasureTextV(WIDER_TEXT, (int)font_size_small);

        text_offset = (int) (1.5 * txt_w_h.y);

        box_w = (int) (txt_w_h.x + 50);

        UpdateCamera(&camera);

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

        handle_input(&draw_config.grid_info.loaded, &ray, &camera);

        if(draw_config.grid_info.loaded ) {

            omp_set_lock(&draw_config.draw_lock);

            if(draw_type == DRAW_FILE) {
                if(draw_config.grid_info.file_name) {
                    window_title = (char *) malloc(strlen(draw_config.grid_info.file_name) + strlen("Visualizing file - ") + 2);
                    sprintf(window_title, "Visualizing file - %s", draw_config.grid_info.file_name);
                    SetWindowTitle(window_title);
                    free(window_title);
                }
            }

            ClearBackground(GRAY);

            BeginMode3D(camera);

            if(!calc_center) {
                if(draw_type == DRAW_SIMULATION) {
                    mesh_offset = find_mesh_center();
                    scale = fmaxf(draw_config.grid_info.grid_to_draw->mesh_side_length.x,
                                  fmaxf(draw_config.grid_info.grid_to_draw->mesh_side_length.y,
                                        draw_config.grid_info.grid_to_draw->mesh_side_length.z)) / 5.0f;
                }
                else {
                    mesh_offset = find_mesh_center_vtk();
                    scale = fmaxf(mesh_offset.x,
                                  fmaxf(mesh_offset.y,
                                        mesh_offset.z)) / 5.0f;
                }
            }

            if(draw_type == DRAW_SIMULATION) {
                draw_alg_mesh(mesh_offset, scale, ray);
            }
            else if(draw_type == DRAW_FILE) {
                draw_vtk_unstructured_grid(mesh_offset, scale, ray);
            }

            omp_unset_lock(&draw_config.draw_lock);

            EndMode3D();

            if(show_scale) {
                draw_scale(font, font_size_small);
            }

            if(hmlen(selected_aps) && show_ap) {
                draw_ap(font, font_size_small, font_size_big);
            }

           if(show_info_box) {
                info_box_h = (text_offset * info_box_lines) + 10;
                draw_box(10, 10, box_w, info_box_h, text_offset, info_box_strings, info_box_lines, font_size_small);
            }

            if(!draw_config.simulating) {
                if(draw_type == DRAW_SIMULATION) {
                    if (show_end_info_box) {
                        configure_end_info_box_strings(&end_info_box_strings);
                        box_h = (text_offset * end_info_box_lines) + 10;
                        draw_box(10, info_box_h + 30, box_w, box_h, text_offset, (const char **) end_info_box_strings,
                                 end_info_box_lines, font_size_small);

                        for (int i = 0; i < end_info_box_lines; i++) {
                            free(end_info_box_strings[i]);
                        }
                    }
                }
            }

            if(show_mesh_info_box) {
                configure_mesh_info_box_strings(&mesh_info_box_strings, draw_type);

                box_h = (text_offset * mesh_info_box_lines) + 10;
                pos_x = GetScreenWidth() - box_w - 10;
                draw_box(pos_x, 10, box_w, box_h, text_offset, (const char **) mesh_info_box_strings,
                         mesh_info_box_lines, font_size_small);

                for (int i = 0; i < mesh_info_box_lines; i++) {
                    free(mesh_info_box_strings[i]);
                }
            }

            if(!draw_config.paused) {
                omp_unset_lock(&draw_config.sleep_lock);
            }

            if(show_selection_box) {
                show_selection_box = !draw_selection_box(font, font_size_small);
            }

            if(show_save_box) {
                show_save_box = !draw_save_box(font, font_size_small);
            }

        }
        else {

            ClearBackground(GRAY);
            float spacing = 20/(float)font.baseSize;
            static Color c = RED;

            if(!draw_config.error_message) {
                draw_config.error_message = strdup("Loading Mesh...");
                c = WHITE;
            }

            error_message_witdh = MeasureTextEx(font, draw_config.error_message, 20, spacing);

            int posx = GetScreenWidth()/2 - (int)error_message_witdh.x/2;
            int posy = GetScreenHeight()/2 - 50;

            int rec_width = (int)(error_message_witdh.x) + 40;

            DrawRectangle(posx, posy, rec_width, 20, c);
            DrawRectangleLines(posx, posy, rec_width, 20, BLACK);
            if(draw_config.error_message) //This should not happen... but it does....
                DrawText(draw_config.error_message, posx + 20, posy, 20, BLACK);

        }

        DrawFPS(GetScreenWidth()  - 100,GetScreenHeight()-20);
        EndDrawing();

    }

    draw_config.exit = true;

    omp_unset_lock(&draw_config.draw_lock);
    omp_unset_lock(&draw_config.sleep_lock);

    CloseWindow();

}
