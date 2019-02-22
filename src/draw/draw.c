//
// Created by sachetto on 11/11/17.
//

#include <float.h>
#include <pthread.h>
#include "draw.h"


#include "../raylib/src/raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "../raylib/src/raygui.h"
#include "../utils/stop_watch.h"
#include "../hash/point_voidp_hash.h"

bool calc_center = false;
//Vector3 *cube_positions = NULL;

int graph_pos_x;
int graph_pos_y;

Font font;

struct ap {
    double v;
    float t;
};

struct point_voidp_hash *selected_aps;

static inline float normalize(float r_min, float r_max, float t_min, float t_max, float m) {
    return ((m - r_min) / (r_max-r_min))*(t_max - t_min) + t_min;
}

static inline Color get_color(double value)
{
    int idx1;        // |-- Our desiBLACK color will be between these two indexes in "color".
    int idx2;        // |
    double fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

    if(value <= 0)      {  idx1 = idx2 = 0;            }    // accounts for an input <=0
    else if(value >= 1)  {  idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
    else
    {
        value = value * (NUM_COLORS-1);        // Will multiply value by 3.
        idx1  = (int)floor(value);                  // Our desiBLACK color will be after this index.
        idx2  = idx1+1;                        // ... and before this index (inclusive).
        fractBetween = value - (double)idx1;    // Distance between the two indexes (0-1).
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

    calc_center = true;

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

    result.x = (max_x+min_x)/2.0f;
    result.y = (max_y+min_y)/2.0f;
    result.z = (max_z+min_z)/2.0f;

    calc_center = true;

    return result;

}

static void draw_alg_mesh(float scale, Ray ray) {

    struct grid *grid_to_draw = draw_config.grid_to_draw;

    Vector3 cubePosition;
    Vector3 cubeSize;
    Color color;

    double max_v = draw_config.max_v;
    double min_v = draw_config.min_v;

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

                struct ap *aps = (struct ap*) point_voidp_hash_search(selected_aps, p);

                if(!draw_config.paused && draw_config.simulating && aps != (void*)-1) {
                    struct ap ap1;
                    ap1.t = draw_config.time;
                    ap1.v = grid_cell->v;
                    sb_push(aps, ap1);
                    point_voidp_hash_insert_or_overwrite(selected_aps, p, aps);
                }

                cubePosition.x = (grid_cell->center_x)/scale;
                cubePosition.y = (grid_cell->center_y)/scale;
                cubePosition.z = (grid_cell->center_z)/scale;


                cubeSize.x = grid_cell->dx/scale;
                cubeSize.y = grid_cell->dy/scale;
                cubeSize.z = grid_cell->dz/scale;

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

                    if(collision) {
                        DrawCubeWiresV(cubePosition, cubeSize, GREEN);

                        if(aps == (void*)-1) {
                            aps = NULL;
                            point_voidp_hash_insert(selected_aps, p, aps);
                            printf("%d\n", selected_aps->n);
                        }
                    }
                }
            }
        }

    }
}

int num_colors = 19;
Color colors[] = {DARKGRAY, YELLOW, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BROWN, DARKBROWN, BLACK, MAGENTA};

void draw_ap() {

    //TODO: draw only when the data is selected
    graph_pos_x = 1;
    graph_pos_y = GetScreenHeight() - 400;

    int grap_width = GetScreenWidth() / 3;

    DrawRectangle(graph_pos_x, graph_pos_y, grap_width, 450, WHITE);
    float min_x = graph_pos_x + 55.0f;
    float max_x = graph_pos_x + grap_width - 10;

    float min_y = graph_pos_y + 350.0f;
    float max_y = graph_pos_y + 50.0f;

    DrawText("Time (ms)", graph_pos_x + 125, min_y + 30, 6, RED);

    Vector2 p1, p2;
    struct ap *aps;
    int c_count = 0;
    for (int i = 0; i < selected_aps->size; i++) {
        for (struct elt *e = selected_aps->table[i % (selected_aps)->size]; e != 0; e = e->next) {
            aps = (struct ap*) e->value;
            int c = sb_count(aps);

            if(c > 0) {
                Color line_color = colors[c_count % num_colors];
                c_count++;
                for (int i = 0; i < c; i++) {

                    if (i + 1 < c) {

                        p1.x = normalize(0.0, draw_config.final_time, min_x, max_x, aps[i].t);
                        p1.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, aps[i].v);

                        p2.x = normalize(0.0, draw_config.final_time, min_x, max_x, aps[i + 1].t);
                        p2.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, aps[i + 1].v);

                        DrawLineV(p1, p2, line_color);
                    }

                }
            }
        }
    }

    int p = 0;
    for(int t = 0; t < draw_config.final_time; ) {
        char tmp[20];


        p1.x = normalize(0.0, draw_config.final_time, min_x, max_x, t);
        p1.y = min_y + 10;

        p2.x = p1.x;
        p2.y = min_y + 15;

        sprintf(tmp, "%d", t);

        if(!(p%2)) {
            DrawTextEx(font, tmp, (Vector2){p1.x - 5, p1.y + 5}, 10, 1, RED);
        }

        p++;

        DrawLineV(p1, p2, RED);

        t += 20;
    }

    int num_ticks = 10;
    float tick_ofsset = (draw_config.max_v - draw_config.min_v)/(float)num_ticks;

    float v = draw_config.min_v;

    for(int t = 0; t <= num_ticks; t++ ) {
        char tmp[20];

        p1.x = graph_pos_x + 35.0f;
        p1.y = normalize(draw_config.min_v, draw_config.max_v, min_y, max_y, v);

        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        sprintf(tmp, "%.2lf", v);
        DrawTextEx(font, tmp, (Vector2){p1.x - 30, p1.y - 5}, 10, 1, RED);
        DrawLineV(p1, p2, RED);

        v += tick_ofsset;
    }

}

void draw_instruction_box () {

    char tmp[100];
    DrawRectangle(10, 10, 320, 193, WHITE);
    DrawRectangleLines(10, 10, 320, 193, BLACK);
    DrawText("Free camera default controls:", 20, 20, 10, BLACK);
    DrawText("- Mouse Wheel to Zoom in-out", 40, 40, 10, DARKGRAY);
    DrawText("- Mouse Wheel Pressed to Pan", 40, 60, 10, DARKGRAY);
    DrawText("- Alt + Mouse Wheel Pressed to Rotate", 40, 80, 10, DARKGRAY);
    DrawText("- Alt + Ctrl + Mouse Wheel Pressed for Smooth Zoom", 40, 100, 10, DARKGRAY);
    DrawText("- Z to reset zoom", 40, 120, 10, DARKGRAY);
    DrawText("- G to only draw the grid lines", 40, 140, 10, DARKGRAY);
    DrawText("- L to enable or disable the grid lines", 40, 160, 10, DARKGRAY);

    sprintf(tmp, "%lf ms", draw_config.time);
    DrawText("Simulation running:", 20, 180, 16, BLACK);
    DrawText(tmp, 170, 180, 16, BLACK);
}

void draw_end_info_box() {

    char tmp[100];

    DrawRectangle(10, 10, 320, 210, SKYBLUE);
    DrawRectangleLines(10, 10, 320, 210, BLUE);

    int text_pos = 20;
    int text_offset = 20;

    DrawText("Simulation finished!", 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "Resolution Time: %ld us\n", draw_config.solver_time);
    DrawText(tmp, 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "ODE Total Time: %ld us\n", draw_config.ode_total_time);
    DrawText(tmp, 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "CG Total Time: %ld us\n", draw_config.cg_total_time);
    DrawText(tmp, 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "Mat time: %ld us\n", draw_config.total_mat_time);
    DrawText(tmp, 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "Refine time: %ld us\n", draw_config.total_ref_time);
    DrawText(tmp, 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "Derefine time: %ld us\n", draw_config.total_deref_time);
    DrawText(tmp, 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "Write time: %ld us\n", draw_config.total_write_time);
    DrawText(tmp, 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "Initial configuration time: %ld us\n", draw_config.total_config_time);
    DrawText(tmp, 20, text_pos, 16, BLACK);

    text_pos += text_offset;

    sprintf(tmp, "CG Total Iterations: %ld\n", draw_config.total_cg_it);
    DrawText(tmp, 20, text_pos, 16, BLACK);
}

const int screenWidth = 1280;
const int screenHeight = 720;

void init_opengl() {

    struct stop_watch timer;

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);

    InitWindow(screenWidth, screenHeight, "Simulation visualization");

    font = LoadFont("misc/Roboto-Black.ttf");

    draw_config.grid_only = false;
    draw_config.grid_lines = true;

    selected_aps = point_voidp_hash_create();

    Camera3D camera;
    camera.position = (Vector3){ 0.064882f, 0.165282f, 15.977825f };  // Camera position
    camera.target = (Vector3){ 0.137565f, 0.199405f, 0.181663f };
    camera.up       = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy     = 45.0f;                                // Camera field-of-view Y
    camera.type     = CAMERA_PERSPECTIVE;                  // Camera mode type

    SetCameraMode(camera, CAMERA_FREE); // Set a free camera mode

    SetTargetFPS(120);

    float scale = 1.0f;

    bool mesh_loaded = false;

    Ray ray = {FLT_MAX, FLT_MAX, FLT_MAX};

    while (!WindowShouldClose()) {

        if (IsKeyDown('Z')) {
            camera.target   = (Vector3){ 0.137565f, 0.199405f, 0.181663f };
        }

        UpdateCamera(&camera);

        if (IsKeyPressed('G')) draw_config.grid_only = !draw_config.grid_only;

        if (IsKeyPressed('L')) draw_config.grid_lines = !draw_config.grid_lines;

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

        if(!draw_config.adaptive) {
            if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                if(!timer.running) {
                    start_stop_watch(&timer);
                }
                else {
                    long delay = stop_stop_watch(&timer);

                    if(delay < 450000) {
                        ray = GetMouseRay(GetMousePosition(), camera);
                    }

                }

            }
        }

        if(draw_config.grid_to_draw && omp_test_lock(&draw_config.draw_lock)) {

            mesh_loaded = true;
            ClearBackground(GRAY);

            BeginMode3D(camera);

            if(!calc_center) {
                scale = fmaxf(draw_config.grid_to_draw->side_length_x, fmaxf(draw_config.grid_to_draw->side_length_y, draw_config.grid_to_draw->side_length_z))/5.0f;
            }

            draw_alg_mesh(scale, ray);
            omp_unset_lock(&draw_config.draw_lock);

            EndMode3D();

            if(draw_config.simulating) {

                draw_ap();
                draw_instruction_box();

            }

            else {
                draw_ap();
                draw_end_info_box();
            }
        }
        else if(!mesh_loaded){
            int posx = GetScreenWidth()/2-150;
            int posy = GetScreenHeight()/2-50;
            ClearBackground(GRAY);
            DrawRectangle(posx, posy, 320, 20, WHITE);
            DrawRectangleLines(posx, posy, 320, 20, BLACK);
            DrawText("Loading Mesh...", posx+80, posy, 20, BLACK);
        }

        if(draw_config.simulating) {
            if(!draw_config.paused) {
                if (GuiButton((Rectangle){ GetScreenWidth() - 170, GetScreenHeight() - 50, 150, 30 }, "Pause Simulation")) {
                    draw_config.paused = true;
                }
            } else {
                if (GuiButton((Rectangle){ GetScreenWidth() - 170, GetScreenHeight() - 50, 150, 30 }, "Continue Simulation")) {
                    draw_config.paused = false;
                }
            }

//            if(!show_ap) {
//                if (GuiButton((Rectangle){ 10, 220, 150, 30 }, "Show AP")) {
//                    show_ap = true;
//                }
//            }
//            else {
//                if (GuiButton((Rectangle){ 10, 220, 150, 30 }, "Hide AP")) {
//                    show_ap = false;
//                }
//            }



            GuiEnable();
        }

        EndDrawing();

    }

    omp_destroy_lock(&draw_config.draw_lock);

    CloseWindow();

}
