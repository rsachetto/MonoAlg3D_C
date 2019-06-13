//
// Created by sachetto on 11/11/17.
//

#include <float.h>
#include <pthread.h>
#include <time.h>


#include "draw.h"
#include "../common_types/common_types.h"
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

static bool show_info_box = true;
static bool show_end_info_box = true;
static bool show_mesh_info_box = true;
static bool show_selection_box = false;

Vector3 max_size;
Vector3 min_size;

#define DOUBLE_CLICK_DELAY 0.5 //seconds
double mouse_timer = -1;

double selected_time = 0.0;
Vector3 current_selected = {FLT_MAX, FLT_MAX, FLT_MAX};

char center_x_text[128] = { 0 };
char center_y_text[128] = { 0 };
char center_z_text[128] = { 0 };

float center_x;
float center_y;
float center_z;


// General variables
Vector2 mousePos = { 0 };
Vector2 windowPos;
bool dragWindow = false;

#define SIZEOF(A) (sizeof(A)/sizeof(A[0]))  // Get number of elements in array `A`. Total size of `A` should be known at compile time.


#define WIDER_TEXT "------------------------------------------------------"

struct point_voidp_hash_entry *selected_aps;

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

    unsigned char red   = (unsigned char) (((color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0]) * 255);
    unsigned char green = (unsigned char) (((color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1]) * 255);
    unsigned char blue  = (unsigned char) (((color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2]) * 255);

    Color result;
    result.r = red;
    result.g = green;
    result.b = blue;
    result.a = 255;

    return result;
}

static bool skip_node(struct cell_node *grid_cell) {

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

    float max_x, max_y, max_z;
    float min_x, min_y, min_z;

    max_x = FLT_MIN;
    max_y = FLT_MIN;
    max_z = FLT_MIN;

    min_x = FLT_MAX;
    min_y = FLT_MAX;
    min_z = FLT_MAX;

    float max_dx, max_dy, max_dz;
    float min_dx, min_dy, min_dz;

    max_dx = FLT_MIN;
    max_dy = FLT_MIN;
    max_dz = FLT_MIN;

    min_dx = FLT_MAX;
    min_dy = FLT_MAX;
    min_dz = FLT_MAX;

    Vector3 result = (Vector3){0.0, 0.0, 0.0};

    if (ac) {
        for (int i = 0; i < n_active; i++) {
            grid_cell = ac[i];
            if(grid_cell->center_x > max_x) {
                max_x = grid_cell->center_x;
                max_dx = grid_cell->dx;
            }
            else if(grid_cell->center_x < min_x) {
                min_x = grid_cell->center_x;
                min_dx = grid_cell->dx;
            }

            if(grid_cell->center_y > max_y) {
                max_y = grid_cell->center_y;
                max_dy = grid_cell->dy;
            }
            else if(grid_cell->center_y < min_y) {
                min_y = grid_cell->center_y;
                min_dy = grid_cell->dy;
            }

            if(grid_cell->center_z > max_z) {
                max_z = grid_cell->center_z;
                max_dz = grid_cell->dz;
            }
            else if(grid_cell->center_z < min_z) {
                min_z = grid_cell->center_z;
                min_dz = grid_cell->dz;
            }

        }
    }

    result.x = max_x;
    if(max_x != min_x)
        result.x = (max_x-min_x)/2.0f;


    result.y = max_y;
    if(max_y != min_y)
        result.y = (max_y-min_y)/2.0f;

    result.z = max_z;
    if(max_z != min_z)
        result.z = (max_z-min_z)/2.0f;

    calc_center = true;

    max_size.x = max_x + max_dx/2.0;
    max_size.y = max_y + max_dy/2.0;
    max_size.z = max_z + max_dz/2.0;

    min_size.x = min_x - min_dx/2.0;
    min_size.y = min_y - min_dy/2.0;
    min_size.z = min_z - min_dz/2.0;

    return result;

}

static Vector3 find_mesh_center_vtk() {

    struct vtk_unstructured_grid *grid_to_draw = draw_config.grid_info.vtk_grid;

    uint32_t n_active = grid_to_draw->num_cells;

    float max_x, max_y, max_z;
    float min_x, min_y, min_z;

    max_x = FLT_MIN;
    max_y = FLT_MIN;
    max_z = FLT_MIN;

    min_x = FLT_MAX;
    min_y = FLT_MAX;
    min_z = FLT_MAX;

    int64_t *cells = grid_to_draw->cells;
    point3d_array points = grid_to_draw->points;

    Vector3 result = (Vector3){0.0, 0.0, 0.0};
    float dx, dy, dz;
    float center_x, center_y, center_z;
    int num_points =grid_to_draw->points_per_cell;

    float max_dx, max_dy, max_dz;
    float min_dx, min_dy, min_dz;

    max_dx = FLT_MIN;
    max_dy = FLT_MIN;
    max_dz = FLT_MIN;

    min_dx = FLT_MAX;
    min_dy = FLT_MAX;
    min_dz = FLT_MAX;

    for (int i = 0; i < n_active*num_points; i+=num_points) {

        dx = fabsf((points[cells[i]].x - points[cells[i+1]].x));
        dy = fabsf((points[cells[i]].y - points[cells[i+3]].y));
        dz = fabsf((points[cells[i]].z - points[cells[i+4]].z));


        center_x = points[cells[i]].x + dx/2.0f;
        center_y = points[cells[i]].y + dy/2.0f;
        center_z = points[cells[i]].z + dz/2.0f;

        if(center_x > max_x) {
            max_x = center_x;
            max_dx = dx;
        }
        else if(center_x < min_x) {
            min_x = center_x;
            min_dx = dx;
        }

        if(center_y > max_y) {
            max_y = center_y;
            max_dy = dy;
        }
        else if(center_y < min_y) {
            min_y = center_y;
            min_dy = dy;
        }

        if(center_z > max_z) {
            max_z = center_z;
            max_dz = dz;
        }
        else if(center_z < min_z) {
            min_z = center_z;
            min_dz = dz;
        }

    }

    result.x = max_x;
    if(max_x != min_x)
        result.x = (max_x-min_x)/2.0f;


    result.y = max_y;
    if(max_y != min_y)
        result.y = (max_y-min_y)/2.0f;

    result.z = max_z;
    if(max_z != min_z)
        result.z = (max_z-min_z)/2.0f;

    calc_center = true;

    max_size.x = max_x + max_dx/2.0;
    max_size.y = max_y + max_dy/2.0;
    max_size.z = max_z + max_dz/2.0;

    min_size.x = min_x - min_dx/2.0;
    min_size.y = min_y - min_dy/2.0;
    min_size.z = min_z - min_dz/2.0;

    return result;

}

Color colors[] = {DARKGRAY, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BROWN, DARKBROWN, BLACK, MAGENTA};
int num_colors = 18;

static void draw_voxel(Vector3 cube_position_draw, Vector3 cube_position_mesh, Vector3 cube_size, real_cpu v, Ray ray) {

    struct point_3d p;
    p.x = cube_position_draw.x;
    p.y = cube_position_draw.y;
    p.z = cube_position_draw.z;

    bool collision;

    real_cpu max_v = draw_config.max_v;
    real_cpu min_v = draw_config.min_v;

    bool grid_only = draw_config.grid_only;
    bool grid_lines = draw_config.grid_lines;
    Color color;

    action_potential_array aps = (struct action_potential*) hmget(selected_aps, p);

    struct action_potential ap1;
    ap1.t = draw_config.time;
    ap1.v = v;

    if(aps != NULL) {
        if(ap1.t > aps[arrlen(aps)-1].t ) {
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

    if(grid_only) {
        DrawCubeWiresV(cube_position_draw, cube_size, color);
    }
    else {

        DrawCubeV(cube_position_draw, cube_size, color);

        if(grid_lines) {
            DrawCubeWiresV(cube_position_draw, cube_size, BLACK);
        }

        if(current_selected.x == cube_position_mesh.x && current_selected.y == cube_position_mesh.y && current_selected.z == cube_position_mesh.z) {
            DrawCubeWiresV(cube_position_draw, cube_size, GREEN);

            if(aps == NULL) {
                arrsetcap(aps, 50);
                struct action_potential ap1;
                ap1.t = draw_config.time;
                ap1.v = v;

                arrput(aps, ap1);
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

    float dx, dy, dz;

    int64_t *cells = grid_to_draw->cells;
    point3d_array points = grid_to_draw->points;

    float center_x, center_y, center_z;
    real_cpu v;

    uint32_t n_active = grid_to_draw->num_cells;

    int num_points = grid_to_draw->points_per_cell;
    int j = num_points;

    for (int i = 0; i < n_active*num_points; i+=num_points) {

        dx = fabsf((points[cells[i]].x - points[cells[i+1]].x));
        dy = fabsf((points[cells[i]].y - points[cells[i+3]].y));
        dz = fabsf((points[cells[i]].z - points[cells[i+4]].z));

        center_x = points[cells[i]].x + dx/2.0f;
        center_y = points[cells[i]].y + dy/2.0f;
        center_z = points[cells[i]].z + dz/2.0f;

        v = grid_to_draw->values[j-num_points];
        j += 1;

        cube_position.x = (float)((center_x - mesh_offset.x)/scale);
        cube_position.y = (float)((center_y - mesh_offset.y)/scale);
        cube_position.z = (float)((center_z - mesh_offset.z)/scale);

        cube_size.x = (float)(dx/scale);
        cube_size.y = (float)(dy/scale);
        cube_size.z = (float)(dz/scale);

        draw_voxel(cube_position, (Vector3){center_x, center_y, center_z}, cube_size, v, ray);

    }
    one_selected = false;
}

static void draw_alg_mesh(Vector3 mesh_offset, real_cpu scale, Ray ray) {

    struct grid *grid_to_draw = draw_config.grid_info.grid_to_draw;

    Vector3 cubePosition;
    Vector3 cubeSize;

    if (grid_to_draw) {

        uint32_t n_active = grid_to_draw->num_active_cells;
        struct cell_node **ac = grid_to_draw->active_cells;
        struct cell_node *grid_cell;

        if (ac) {
            for (int i = 0; i < n_active; i++) {

                grid_cell = ac[i];

                if(skip_node(grid_cell)) {
                    continue;
                }


                cubePosition.x = (float)((grid_cell->center_x - mesh_offset.x)/scale);
                cubePosition.y = (float)((grid_cell->center_y - mesh_offset.y)/scale);
                cubePosition.z = (float)((grid_cell->center_z - mesh_offset.z)/scale);

                cubeSize.x = (float)(grid_cell->dx/scale);
                cubeSize.y = (float)(grid_cell->dy/scale);
                cubeSize.z = (float)(grid_cell->dz/scale);

                draw_voxel(cubePosition, (Vector3){grid_cell->center_x, grid_cell->center_y, grid_cell->center_z}, cubeSize, ac[i]->v, ray);

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

void draw_ap(Font font, int font_size_small, int font_size_big) {

    int graph_height = GetScreenHeight()/3;

    int graph_pos_x = 1;
    int graph_pos_y = GetScreenHeight() - graph_height;

    int graph_width = GetScreenWidth();

    float spacing_big = font_size_big/(float)font.baseSize;
    float spacing_small = font_size_small/(float)font.baseSize;

    DrawRectangle(graph_pos_x, graph_pos_y, graph_width, graph_height, WHITE);

    char tmp[256];
    sprintf(tmp, "%.2lf", draw_config.final_time);

    Vector2 width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    float min_x = graph_pos_x + 2.0f*width.x;
    float max_x = graph_pos_x + graph_width - width.x;

    float min_y = graph_pos_y + graph_height - 30.0f;
    float max_y = graph_pos_y + 20.0f; //This is actually the smallest allowed y

    int n = hmlen(selected_aps);

    char *ap_text = "%d AP(s) selected ( cell at %f, %f, %f )";

    if(draw_selected_ap_text) {
        double time_elapsed = GetTime() - selected_time;
        unsigned char alpha = (unsigned char) clamp(255 - time_elapsed*25, 0, 255);
        Color c = colors[(n-1) % num_colors];
        c.a = alpha;

        sprintf(tmp, ap_text, n, current_selected.x, current_selected.y, current_selected.z);

        width = MeasureTextEx(font, ap_text, font_size_big, spacing_big);

        DrawTextEx(font, tmp, (Vector2){graph_pos_x + graph_width/2.0f - width.x/2.0f - min_x, max_y}, font_size_big, 1, c);

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

    DrawTextEx(font, time_text, (Vector2){min_x + graph_width/2.0f - width.x/2.0f, (float)min_y + 50.0f}, font_size_big, spacing_big, BLACK);

    Vector2 p1, p2;


    real_cpu time = 0.0;

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
    Vector2 max_w = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    float graph_min_y = min_y;

    //Draw vertical ticks
    for(int t = 0; t <= num_ticks; t++ ) {

        p1.x = graph_pos_x + 5.0f;
        p1.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, v);

        if(t == 0) {
            graph_min_y = p1.y;
        }

        sprintf(tmp, "%.2lf", v);
        width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

        DrawTextEx(font, tmp, (Vector2){p1.x + (max_w.x - width.x), p1.y - width.y/2.0f}, font_size_small, spacing_small, RED);

        p1.x = p1.x + max_w.x + 5.0f;
        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        DrawLineV(p1, (Vector2){max_x, p1.y}, LIGHTGRAY);

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


    float graph_min_x = ((graph_pos_x + 5.0f + max_w.x + 5.0f) + (graph_pos_x + 5.0f + max_w.x + 5.0f + 10.0))/2.0;

    //Draw horizontal ticks
    for(int t = 0; t <= num_ticks; t++) {

        if(t == 0) {
            p1.x =  graph_min_x;
        }
        else {
            p1.x = normalize(0.0, draw_config.final_time, min_x, max_x, time) - 20;

        }
        p1.y = graph_min_y - 5;

        p2.x = p1.x;
        p2.y = graph_min_y + 5;

        if(!(t%2)) {
            sprintf(tmp, "%.2lf", time);
            width = MeasureTextEx(font, tmp, font_size_small, spacing_small);
            DrawTextEx(font, tmp, (Vector2){p1.x - width.x/2.0f, p1.y + 10}, font_size_small, spacing_small, RED);
        }

        DrawLineV(p1, p2, RED);

        DrawLineV(p1, (Vector2){p1.x, max_y}, LIGHTGRAY);


        time+= tick_ofsset;
    }

    //Draw vertical line
    p1.x = ((graph_pos_x + 5.0f + max_w.x + 5.0f) + (graph_pos_x + 5.0f + max_w.x + 5.0f + 10.0))/2.0 ;
    p1.y = min_y + 5;

    p2.x = p1.x;
    p2.y  = max_y;
    DrawLineV(p1, p2, RED);

    //Draw horizontal line
    p1.x = min_x - 20;
    p1.y = graph_min_y;

    p2.x = max_x;
    p2.y = p1.y;
    DrawLineV(p1, p2, RED);


    struct action_potential *aps;

    for (int j = 0; j < n; j++) {

        aps = (struct action_potential*) selected_aps[j].value;
        int c = arrlen(aps);

        if(c > 0) {
            Color line_color = colors[j % num_colors];
            for (int i = 0; i < c; i++) {

                if(aps[i].t <= draw_config.time) {
                    if (i + 1 < c) {

                        p1.x = normalize(0.0, draw_config.final_time, graph_min_x, max_x, aps[i].t);
                        p1.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, aps[i].v);

                        p2.x = normalize(0.0, draw_config.final_time, graph_min_x, max_x, aps[i + 1].t);
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
}

static void draw_scale(Font font, int font_size_small) {

    float initial_y =  GetScreenHeight()/2.0;

    float spacing_small = font_size_small/(float)font.baseSize;

    int num_ticks;
    real_cpu tick_ofsset = 12;
    num_ticks = (int) ((draw_config.max_v - draw_config.min_v)/tick_ofsset);

    if(num_ticks < 4) {
        num_ticks = 4;
        tick_ofsset = (draw_config.max_v - draw_config.min_v)/num_ticks;
    }

    char tmp[256];

    real_cpu v = draw_config.min_v;
    sprintf(tmp, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    Vector2 p1, p2, width;
    float scale_pos_x = GetScreenWidth()-30.0;

    float scale_rec_size = 30.0;
    Color color;

    real_cpu  min_v = draw_config.min_v;
    real_cpu  max_v = draw_config.max_v;

    for(int t = 0; t <= num_ticks; t++ ) {

        p1.x = scale_pos_x - 55.0f;
        p1.y = initial_y + scale_rec_size/2.0;

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

    sprintf(tmp, "Resolution Time: %lf s", draw_config.solver_time/1000.0/1000.0);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "ODE Total Time: %lf s", draw_config.ode_total_time/1000.0/1000.0);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "CG Total Time: %lf s", draw_config.cg_total_time/1000.0/1000.0);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Mat time: %lf s", draw_config.total_mat_time/1000.0/1000.0);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Refine time: %lf s", draw_config.total_ref_time/1000.0/1000.0);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Derefine time: %lf s", draw_config.total_deref_time/1000.0/1000.0);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Write time: %lf s", draw_config.total_write_time/1000.0/1000.0);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, "Initial configuration time: %lf s", draw_config.total_config_time/1000.0/1000.0);
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
//    float sx, sy, sz;

    if(draw_type == DRAW_SIMULATION) {
        n_active = draw_config.grid_info.grid_to_draw->num_active_cells;
    }
    else {
        n_active = draw_config.grid_info.vtk_grid->num_cells;
    }

    (*(info_string))[index++] = strdup("Mesh information:");

    sprintf(tmp, " - Num. of Volumes: %d", n_active);
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

    int pos_x = (int) windowPos.x;
    int pos_y = (int) windowPos.y;

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

static void handle_input(bool *mesh_loaded, Ray *ray, Camera3D *camera) {

    {
        mousePos = GetMousePosition();

        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON))
        {
            if (CheckCollisionPointRec(mousePos, (Rectangle){  windowPos.x,  windowPos.y, box_width - 18 , WINDOW_STATUSBAR_HEIGHT }))
            {
                dragWindow = true;
            }
        }

        if (dragWindow)
        {
            windowPos.x = (mousePos.x) - (box_width - 18)/2;
            windowPos.y = (mousePos.y) - WINDOW_STATUSBAR_HEIGHT/2;

            if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) dragWindow = false;

        }
    }

    if(IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown((KEY_LEFT_CONTROL))) {
        if(IsKeyPressed(KEY_F)) {
            show_selection_box = true;
            show_selection_box = true;
            windowPos.x = GetScreenWidth() / 2 - box_width;
            windowPos.y = GetScreenHeight() / 2 - box_height;
        }
    }

    if (IsKeyPressed('S')) {
        show_scale = !show_scale;
        return;
    }

    if(draw_config.paused) {

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
        draw_config.grid_only = !draw_config.grid_only;
        return;
    }

    if (IsKeyPressed('L')) {
        draw_config.grid_lines = !draw_config.grid_lines;
        return;
    }

    if (IsKeyPressed(KEY_SPACE)) {
        draw_config.paused = !draw_config.paused;
        return;
    }

    if (IsKeyPressed('R')) {
        for(int i = 0; i < hmlen(selected_aps); i++) {
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

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    char *tmp = NULL;

    int draw_type = draw_config.draw_type;

    if(draw_type == DRAW_SIMULATION) {
        tmp = (char *) malloc(strlen(draw_config.config_name) + strlen("Simulation visualization - ") + 2);
        sprintf(tmp, "Simulation visualization - %s", draw_config.config_name);
    }
    else {
        tmp = strdup("Opening mesh...");
    }

    InitWindow(0, 0, tmp);

    free(tmp);

    int current_monitor = 0;

    Font font = GetFontDefault();

    bool mesh_loaded = false;

    draw_config.grid_only = false;
    draw_config.grid_lines = true;

    selected_aps = NULL;
    hmdefault(selected_aps, NULL);

    Camera3D camera;

    camera.position = (Vector3){  11.082402f, 6.763101, 8.921088  };  // Camera position
    camera.target = (Vector3){ 1.081274, -1.581945f, 0.196326};
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
            " - S to show/hide scale",
            " - C to show/hide everything except grid",
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
    int font_size_small;
    int font_size_big;

    end_info_box_strings = (char **) malloc(sizeof(char *) * end_info_box_lines);
    mesh_info_box_strings = (char **) malloc(sizeof(char *) * mesh_info_box_lines);
    bool have_grid;

    Vector2 error_message_witdh;

    while (!WindowShouldClose()) {
        //Configure font size according to monitor resolution
        current_monitor = GetCurrentMonitor();

        font_size_small  = (int)(10*(GetMonitorHeight(current_monitor)/1080.0));
        if(font_size_small < 10) font_size_small = 10;

        font_size_big    = (int)(16*(GetMonitorHeight(current_monitor)/1080.0));
        if(font_size_big < 16) font_size_big = 16;

        txt_w_h = MeasureTextV(WIDER_TEXT, font_size_small);

        text_offset = (int) (1.5 * txt_w_h.y);

        box_w = (int) (txt_w_h.x + 50);

        UpdateCamera(&camera);

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

        handle_input(&mesh_loaded, &ray, &camera);

        if(draw_type == DRAW_SIMULATION) {
            have_grid = draw_config.grid_info.grid_to_draw;
        } else {
            have_grid = draw_config.grid_info.vtk_grid;
        }

        if(have_grid) {

            omp_set_lock(&draw_config.draw_lock);

            if(draw_type == DRAW_FILE) {
                if(draw_config.grid_info.file_name) {
                    tmp = (char *) malloc(strlen(draw_config.grid_info.file_name) + strlen("Visualizing file - ") + 2);
                    sprintf(tmp, "Visualizing file - %s", draw_config.grid_info.file_name);
                    SetWindowTitle(tmp);
                    free(tmp);
                }
            }

            mesh_loaded = true;
            ClearBackground(GRAY);

            BeginMode3D(camera);

            if(!calc_center) {
                if(draw_type == DRAW_SIMULATION) {
                    mesh_offset = find_mesh_center();
                    scale = fmaxf(draw_config.grid_info.grid_to_draw->side_length_x,
                                  fmaxf(draw_config.grid_info.grid_to_draw->side_length_y,
                                        draw_config.grid_info.grid_to_draw->side_length_z)) / 5.0f;
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


        }
        else if(!mesh_loaded) {

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
