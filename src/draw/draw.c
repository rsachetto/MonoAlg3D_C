//
// Created by sachetto on 11/11/17.
//

#include <float.h>
#include <pthread.h>
#include "draw.h"


#include "../raylib/src/raylib.h"

#include "../common_types/common_types.h"

#include "../single_file_libraries/stb_ds.h"
#include "../raylib/src/camera.h"

static bool calc_center = false;
static bool one_selected = false;
static bool draw_selected_ap_text = false;
static bool show_ap = true;
static bool show_scale = true;
static bool c_pressed = false;

static bool show_info_box = true;
static bool show_end_info_box = true;
static bool show_mesh_info_box = true;

#define DOUBLE_CLICK_DELAY 0.5 //seconds

double selected_time = 0.0;

#define WIDER_TEXT "------------------------------------------------------"

struct action_potential {
    real_cpu v;
    real_cpu t;
};

typedef struct action_potential * action_potential_array;

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

static bool have_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell) {

    struct cell_node *black_neighbor_cell;

    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    struct transition_node *white_neighbor_cell;

    if(neighbour_grid_cell_level > grid_cell->cell_data.level)
    {
        if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE)
        {
            while(true)
            {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE)
                {
                    white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                    if(white_neighbor_cell->single_connector == NULL)
                    {
                        return false;
                    }
                    else
                    {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                }
                else
                {
                    black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
                    return black_neighbor_cell->active;
                }
            }
        }
        else {
            black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
            return black_neighbor_cell->active;
        }

    }
    else
    {
        if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE)
        {
            while(true)
            {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE)
                {
                    white_neighbor_cell = (struct transition_node *)(neighbour_grid_cell);
                    if(white_neighbor_cell->single_connector == NULL)
                    {
                        return false;
                    }
                    else
                    {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                }
                else
                {
                    black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
                    return black_neighbor_cell->active;
                }
            }
        }
        else {
            black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
            return black_neighbor_cell->active;
        }
    }
}

static bool skip_node(struct cell_node *grid_cell) {

    if(!have_neighbour(grid_cell, grid_cell->north) ) {
        return false;
    }
    else if(!have_neighbour(grid_cell, grid_cell->south) ) {
        return false;
    }
    else if(!have_neighbour(grid_cell, grid_cell->west) ) {
        return false;
    }
    else if(!have_neighbour(grid_cell, grid_cell->east) ) {
        return false;
    }
    else if(!have_neighbour(grid_cell, grid_cell->front) ) {
        return false;
    }
    else if(!have_neighbour(grid_cell, grid_cell->back) ) {
        return false;
    }
    else {
        return true;
    }
}

static Vector3 find_mesh_center() {

    struct grid *grid_to_draw = draw_config.grid_to_draw;

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

    Vector3 result = (Vector3){0.0, 0.0, 0.0};

    if (ac) {
        for (int i = 0; i < n_active; i++) {
            grid_cell = ac[i];
            if(grid_cell->center_x > max_x) {
                max_x = grid_cell->center_x;
            }
            else if(grid_cell->center_x < min_x) {
                min_x = grid_cell->center_x;
            }

            if(grid_cell->center_y > max_y) {
                max_y = grid_cell->center_y;
            }
            else if(grid_cell->center_y < min_y) {
                min_y = grid_cell->center_y;
            }

            if(grid_cell->center_z > max_z) {
                max_z = grid_cell->center_z;
            }
            else if(grid_cell->center_z < min_z) {
                min_z = grid_cell->center_z;
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

    return result;

}

static Vector3 find_mesh_center_vtk() {

    struct vtk_unstructured_grid *grid_to_draw = draw_config.vtk_grid;

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

    for (int i = 0; i < n_active*num_points; i+=num_points) {


        dx = fabsf((points[cells[i]].x - points[cells[i+1]].x));
        dy = fabsf((points[cells[i]].y - points[cells[i+3]].y));
        dz = fabsf((points[cells[i]].z - points[cells[i+4]].z));


        center_x = points[cells[i]].x + dx/2.0f;
        center_y = points[cells[i]].y + dy/2.0f;
        center_z = points[cells[i]].z + dz/2.0f;

        if(center_x > max_x) {
            max_x = center_x;
        }
        else if(center_x < min_x) {
            min_x = center_x;
        }

        if(center_y > max_y) {
            max_y = center_y;
        }
        else if(center_y < min_y) {
            min_y = center_y;
        }

        if(center_z > max_z) {
            max_z = center_z;
        }
        else if(center_z < min_z) {
            min_z = center_z;
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

    return result;

}

Color colors[] = {DARKGRAY, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BROWN, DARKBROWN, BLACK, MAGENTA};
int num_colors = 18;

static void draw_vtk_unstructured_grid(Vector3 mesh_offset, real_cpu scale, Ray ray) {

    struct vtk_unstructured_grid *grid_to_draw = draw_config.vtk_grid;

    Vector3 cubePosition;
    Vector3 cubeSize;
    Color color;

    real_cpu max_v = draw_config.max_v;
    real_cpu min_v = draw_config.min_v;

    bool grid_only = draw_config.grid_only;
    bool grid_lines = draw_config.grid_lines;

    bool collision;
    int64_t *cells = grid_to_draw->cells;
    point3d_array points = grid_to_draw->points;

    float dx, dy, dz;
    float center_x, center_y, center_z;
    real_cpu v;



    if (grid_to_draw) {

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

//            if(skip_node(grid_cell)) {
//                continue;
//            }

            struct point_3d p;
            p.x = center_x;
            p.y = center_y;
            p.z = center_z;

            action_potential_array aps = (struct action_potential*) hmget(selected_aps, p);

            if(!draw_config.paused && draw_config.simulating && aps != NULL) {
                struct action_potential ap1;
                ap1.t = draw_config.time;
                ap1.v = v;
                arrput(aps, ap1);

                hmput(selected_aps, p, aps);
            }

            cubePosition.x = (float)((center_x - mesh_offset.x)/scale);
            cubePosition.y = (float)((center_y - mesh_offset.y)/scale);
            cubePosition.z = (float)((center_z - mesh_offset.z)/scale);

            cubeSize.x = (float)(dx/scale);
            cubeSize.y = (float)(dy/scale);
            cubeSize.z = (float)(dz/scale);

            collision = CheckCollisionRayBox(ray,
                                             (BoundingBox){(Vector3){ cubePosition.x - cubeSize.x/2, cubePosition.y - cubeSize.y/2, cubePosition.z - cubeSize.z/2 },
                                                           (Vector3){ cubePosition.x + cubeSize.x/2, cubePosition.y + cubeSize.y/2, cubePosition.z + cubeSize.z/2 }});

            color = get_color((v - min_v)/(max_v - min_v));

            if(grid_only) {
                DrawCubeWiresV(cubePosition, cubeSize, color);
            }
            else {

                DrawCubeV(cubePosition, cubeSize, color);

                if(grid_lines) {
                    DrawCubeWiresV(cubePosition, cubeSize, BLACK);
                }

                if(collision && !one_selected) {
                    DrawCubeWiresV(cubePosition, cubeSize, GREEN);

                    if(aps == NULL) {
                        arrsetcap(aps, 50);
                        hmput(selected_aps, p, aps);
                        draw_selected_ap_text = true;
                        selected_time = GetTime();
                    }
                    one_selected = true;

                }
            }
        }
    }
    one_selected = false;
}

static void draw_alg_mesh(Vector3 mesh_offset, real_cpu scale, Ray ray) {

    struct grid *grid_to_draw = draw_config.grid_to_draw;

    Vector3 cubePosition;
    Vector3 cubeSize;
    Color color;

    real_cpu max_v = draw_config.max_v;
    real_cpu min_v = draw_config.min_v;

    bool grid_only = draw_config.grid_only;
    bool grid_lines = draw_config.grid_lines;

    bool collision;

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

                struct point_3d p;
                p.x = grid_cell->center_x;
                p.y = grid_cell->center_y;
                p.z = grid_cell->center_z;

                action_potential_array aps = (struct action_potential*) hmget(selected_aps, p);

                if(!draw_config.paused && draw_config.simulating && aps != NULL) {
                    struct action_potential ap1;
                    ap1.t = draw_config.time;
                    ap1.v = grid_cell->v;
                    arrput(aps, ap1);

                    hmput(selected_aps, p, aps);
                }

                cubePosition.x = (float)((grid_cell->center_x - mesh_offset.x)/scale);
                cubePosition.y = (float)((grid_cell->center_y - mesh_offset.y)/scale);
                cubePosition.z = (float)((grid_cell->center_z - mesh_offset.z)/scale);

                cubeSize.x = (float)(grid_cell->dx/scale);
                cubeSize.y = (float)(grid_cell->dy/scale);
                cubeSize.z = (float)(grid_cell->dz/scale);

                collision = CheckCollisionRayBox(ray,
                                                 (BoundingBox){(Vector3){ cubePosition.x - cubeSize.x/2, cubePosition.y - cubeSize.y/2, cubePosition.z - cubeSize.z/2 },
                                                               (Vector3){ cubePosition.x + cubeSize.x/2, cubePosition.y + cubeSize.y/2, cubePosition.z + cubeSize.z/2 }});

                color = get_color((grid_cell->v - min_v)/(max_v - min_v));

                if(grid_only) {
                    DrawCubeWiresV(cubePosition, cubeSize, color);
                }
                else {

                    DrawCubeV(cubePosition, cubeSize, color);

                    if(grid_lines) {
                        DrawCubeWiresV(cubePosition, cubeSize, BLACK);
                    }

                    if(collision && !one_selected) {
                        DrawCubeWiresV(cubePosition, cubeSize, GREEN);

                        if(aps == NULL) {
                            arrsetcap(aps, 50);
                            hmput(selected_aps, p, aps);
                            draw_selected_ap_text = true;
                            selected_time = GetTime();
                        }
                        one_selected = true;

                    }
                }
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

    char *ap_text = "%d AP(s) selected";

    width = MeasureTextEx(font, ap_text, font_size_big, spacing_big);

    if(draw_selected_ap_text) {
        double time_elapsed = GetTime() - selected_time;
        unsigned char alpha = (unsigned char) clamp(255 - time_elapsed*100, 0, 255);
        Color c = colors[(n-1) % num_colors];
        c.a = alpha;
        sprintf(tmp, ap_text, n);
        DrawTextEx(font, tmp, (Vector2){graph_pos_x + graph_width/2.0f - width.x/2.0f - min_x, max_y}, font_size_big, 1, c);

        if(alpha == 0) {
            draw_selected_ap_text = false;
            selected_time = 0.0;
        }
    }

    char *time_text = "Time (ms)";
    width = MeasureTextEx(font, time_text, font_size_big, spacing_big);

    DrawTextEx(font, time_text, (Vector2){min_x + graph_width/2.0f - width.x/2.0f, (float)min_y + 50.0f}, font_size_big, spacing_big, BLACK);

    Vector2 p1, p2;


    real_cpu time = 0.0;

    int num_ticks;
    real_cpu tick_ofsset = 10;
    num_ticks = (draw_config.max_v - draw_config.min_v)/tick_ofsset;

    if(num_ticks < 4) {
        num_ticks = 4;
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
    num_ticks = draw_config.final_time/tick_ofsset;

    if(num_ticks < 4) {
        num_ticks = 4;
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

                if (i + 1 < c) {

                    p1.x = normalize(0.0, draw_config.final_time, graph_min_x, max_x, aps[i].t);
                    p1.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, aps[i].v);

                    p2.x = normalize(0.0, draw_config.final_time, graph_min_x, max_x, aps[i + 1].t);
                    p2.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, aps[i + 1].v);

                    //TODO: create an option for this???
                    if(aps[i+1].v > draw_config.max_v) draw_config.max_v = aps[i+1].v;
                    if(aps[i+1].v < draw_config.min_v) draw_config.min_v = aps[i+1].v;

                    DrawLineV(p1, p2, line_color);
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
    num_ticks = (draw_config.max_v - draw_config.min_v)/tick_ofsset;

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


        DrawRectangle(scale_pos_x,initial_y, 20, 30, color);
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

    (*(info_string))[index++] = strdup("Simulation finished!");

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
    float sx, sy, sz;

    if(draw_type == DRAW_SIMULATION) {
        n_active = draw_config.grid_to_draw->num_active_cells;
        sx = draw_config.grid_to_draw->side_length_x;
        sy = draw_config.grid_to_draw->side_length_y;
        sz = draw_config.grid_to_draw->side_length_z;
    }
    else {
        n_active = draw_config.vtk_grid->num_cells;
        sx = 0.0;
        sy = 0.0;
        sz = 0.0;
    }

    (*(info_string))[index++] = strdup("Mesh information:");

    sprintf(tmp, " - Num. of Volumes: %d", n_active);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Max X: %f", sx);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Max Y: %f", sy);
    (*(info_string))[index++]  = strdup(tmp);

    sprintf(tmp, " - Max Z: %f", sz);
    (*(info_string))[index++]  = strdup(tmp);

    if(draw_config.paused) {
        sprintf(tmp, "Simulation paused: %lf of %lf ms", draw_config.time, draw_config.final_time);
    }
    else if(draw_config.simulating){
        sprintf(tmp, "Simulation running: %lf of %lf ms", draw_config.time, draw_config.final_time);

    }
    else  {
        sprintf(tmp, "Simulation finished: %lf of %lf ms", draw_config.time, draw_config.final_time);
    }

    (*(info_string))[index]  = strdup(tmp);


}

void handle_input(bool *mesh_loaded, Ray *ray) {
    if (IsKeyPressed('G')) draw_config.grid_only = !draw_config.grid_only;

    if (IsKeyPressed('L')) draw_config.grid_lines = !draw_config.grid_lines;

    if (IsKeyPressed(KEY_SPACE)) {

        if(draw_config.paused) {
            omp_unset_lock(&draw_config.sleep_lock);
        }

        draw_config.paused = !draw_config.paused;
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
        draw_config.grid_to_draw = NULL;
        *mesh_loaded = false;

        ray->position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
        ray->direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
        calc_center = false;
    }

    if (IsKeyPressed('A')) {
        show_ap = !show_ap;
    }


    if (IsKeyPressed('S')) {
        show_scale = !show_scale;
    }

    if (IsKeyPressed('C')) {

        show_scale = c_pressed;
        show_ap = c_pressed;
        show_info_box = c_pressed;
        show_end_info_box = c_pressed;
        show_mesh_info_box = c_pressed;

        c_pressed = !c_pressed;
    }


}

void init_and_open_visualization_window(int draw_type) {

    omp_set_lock(&draw_config.sleep_lock);

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    char *tmp = NULL;

    if(draw_type == DRAW_SIMULATION) {
        tmp = (char *) malloc(strlen(draw_config.config_name) + strlen("Simulation visualization - ") + 2);
        sprintf(tmp, "Simulation visualization - %s", draw_config.config_name);
    }
    else if(draw_type == DRAW_FILE) {
        tmp = strdup("Visualizing...");
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

    double mouse_timer = -1;

    char **end_info_box_strings = NULL;
    char **mesh_info_box_strings = NULL;

    const char *info_box_strings[] = {
            "Default controls:",
            " - Mouse Wheel to Zoom in-out",
            " - Mouse Wheel Pressed to Pan",
            " - Alt + Mouse Wheel Pressed to Rotate",
            " - Alt + Ctrl + Mouse Wheel Pressed for Smooth Zoom",
            " - Z to reset zoom",
            " - G to only draw the grid lines",
            " - L to enable or disable the grid lines",
            " - Double click on a volume to show the AP",
            " - R to restart simulation",
            " - A to show/hide AP visualization",
            " - S to show/hide scale",
            " - C to show/hide everything except grid",
            " - Space to start or pause simulation"
    };

    int info_box_lines = sizeof(info_box_strings)/sizeof(info_box_strings[0]);
    int end_info_box_lines = 11;
    int mesh_info_box_lines = 6;

    Vector2 txt_w_h;

    int text_offset;

    int box_w, box_h;

    int pos_x;
    int font_size_small;
    int font_size_big;

    end_info_box_strings = (char **) malloc(sizeof(char *) * end_info_box_lines);
    mesh_info_box_strings = (char **) malloc(sizeof(char *) * mesh_info_box_lines);
    bool have_grid = false;

    while (!WindowShouldClose()) {

        current_monitor = GetCurrentMonitor();

        font_size_small  = (int)(10*(GetMonitorHeight(current_monitor)/1080.0));
        if(font_size_small < 10) font_size_small = 10;

        font_size_big    = (int)(16*(GetMonitorHeight(current_monitor)/1080.0));
        if(font_size_big < 10) font_size_big = 16;

        txt_w_h = MeasureTextV(WIDER_TEXT, font_size_small);

        text_offset = (int) (1.5 * txt_w_h.y);

        box_w = (int) (txt_w_h.x + 50);

        if (IsKeyDown('Z')) {
            camera.target = (Vector3){ 1.081274, -1.581945f, 0.196326};;
        }

        UpdateCamera(&camera);

        handle_input(&mesh_loaded, &ray);

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            if(mouse_timer == -1) {
                mouse_timer = GetTime();
            }
            else {
                double delay = GetTime()-mouse_timer;

                if(delay < DOUBLE_CLICK_DELAY) {
                    ray = GetMouseRay(GetMousePosition(), camera);
                    mouse_timer = -1;
                }
                else {
                    mouse_timer = -1;
                }

            }
        }

        if(draw_type == DRAW_SIMULATION) {
            have_grid = draw_config.grid_to_draw;
        } else {
            have_grid = draw_config.vtk_grid;
        }

        if(have_grid) {

            omp_set_lock(&draw_config.draw_lock);

            mesh_loaded = true;
            ClearBackground(GRAY);

            BeginMode3D(camera);

            if(!calc_center) {
                if(draw_type == DRAW_SIMULATION) {
                    mesh_offset = find_mesh_center();
                    scale = fmaxf(draw_config.grid_to_draw->side_length_x,
                                  fmaxf(draw_config.grid_to_draw->side_length_y,
                                        draw_config.grid_to_draw->side_length_z)) / 5.0f;
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

            if(draw_config.simulating) {
                if(show_info_box) {
                    box_h = (text_offset * info_box_lines) + 10;
                    draw_box(10, 10, box_w, box_h, text_offset, info_box_strings, info_box_lines, font_size_small);
                }
            }

            else {
                if(show_end_info_box) {
                    configure_end_info_box_strings(&end_info_box_strings);
                    box_h = (text_offset * end_info_box_lines) + 10;
                    draw_box(10, 10, box_w, box_h, text_offset, (const char **) end_info_box_strings,
                             end_info_box_lines, font_size_small);

                    for (int i = 0; i < end_info_box_lines; i++) {
                        free(end_info_box_strings[i]);
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

        }
        else if(!mesh_loaded) {
            int posx = GetScreenWidth()/2 - 150;
            int posy = GetScreenHeight()/2 - 50;
            ClearBackground(GRAY);
            DrawRectangle(posx, posy, 320, 20, WHITE);
            DrawRectangleLines(posx, posy, 320, 20, BLACK);
            DrawText("Loading Mesh...", posx+80, posy, 20, BLACK);
        }

        EndDrawing();

    }

    omp_destroy_lock(&draw_config.draw_lock);
    omp_destroy_lock(&draw_config.sleep_lock);

    if(draw_config.paused) {
        omp_unset_lock(&draw_config.sleep_lock);
        draw_config.paused = false;
    }

    draw_config.exit = true;
    CloseWindow();

}
