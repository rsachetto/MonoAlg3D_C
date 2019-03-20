//
// Created by sachetto on 11/11/17.
//

#include <float.h>
#include <pthread.h>
#include "draw.h"


#include "../raylib/src/raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "../raylib/src/raygui.h"
#include "../common_types/common_types.h"

#include "../single_file_libraries/stb_ds.h"

static bool calc_center = false;
static bool one_selected = false;
static bool draw_selected_ap_text = false;

#define DOUBLE_CLICK_DELAY 0.5 //seconds

double selected_time = 0.0;

int graph_pos_x;
int graph_pos_y;

Font font;

struct action_potential {
    real_cpu v;
    real_cpu t;
};

typedef struct action_potential * action_potential_array;

struct point_voidp_hash_entry *selected_aps;

static inline real_cpu normalize(real_cpu r_min, real_cpu r_max, real_cpu t_min, real_cpu t_max, real_cpu m) {
    return ((m - r_min) / (r_max-r_min))*(t_max - t_min) + t_min;
}

static inline Color get_color(real_cpu value)
{
    int idx1;        // |-- Our desiBLACK color will be between these two indexes in "color".
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

    real_cpu max_x, max_y, max_z;
    real_cpu min_x, min_y, min_z;

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

    result.x = (float)(max_x+min_x)/2.0f;
    result.y = (float)(max_y+min_y)/2.0f;
    result.z = (float)(max_z+min_z)/2.0f;

    calc_center = true;

    return result;

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

int num_colors = 19;
Color colors[] = {DARKGRAY, YELLOW, GOLD, ORANGE, PINK, RED, MAROON, GREEN, LIME, DARKGREEN, BLUE, DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BROWN, DARKBROWN, BLACK, MAGENTA};

 double clamp(double x, double min, double max) {
     if (x < min)
         x = min;
     else if (x > max)
         x = max;
     return x;
 }

void draw_ap() {

    graph_pos_x = 1;
    graph_pos_y = GetScreenHeight() - 400;

    int grap_width = GetScreenWidth() / 3;

    DrawRectangle(graph_pos_x, graph_pos_y, grap_width, 450, WHITE);
    real_cpu min_x = graph_pos_x + 55.0f;
    real_cpu max_x = graph_pos_x + grap_width - 10;

    real_cpu min_y = graph_pos_y + 350.0f;
    real_cpu max_y = graph_pos_y + 50.0f;

    int n = hmlen(selected_aps);

    if(draw_selected_ap_text) {
        char tmp[256];
        double time_elapsed = GetTime() - selected_time;
        unsigned char alpha = (unsigned char) clamp(255 - time_elapsed*100, 0, 255);
        Color c = colors[(n-1) % num_colors];
        c.a = alpha;
        sprintf(tmp, "%d AP(s) selected", n);
        DrawTextEx(font, tmp, (Vector2){graph_pos_x + 160.0f, (float)max_y}, 16, 1, c);

        if(alpha == 0) {
            draw_selected_ap_text = false;
            selected_time = 0.0;
        }
    }

    DrawTextEx(font, "Time (ms)", (Vector2){graph_pos_x + 160.0f, (float)min_y + 30.0f}, 16, 1, BLACK);

    Vector2 p1, p2;
    struct action_potential *aps;

    for (int j = 0; j < n; j++) {

        aps = (struct action_potential*) selected_aps[j].value;
        int c = arrlen(aps);

        if(c > 0) {
            Color line_color = colors[j % num_colors];
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
    real_cpu tick_ofsset = (draw_config.max_v - draw_config.min_v)/(real_cpu)num_ticks;

    real_cpu v = draw_config.min_v;

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

    int text_position = 10;
    int text_offset = 20;

    int box_w = 320;
    int box_h = 263;

    char tmp[100];
    DrawRectangle(10, 10, box_w, box_h, WHITE);

    DrawRectangleLines(10, 10, box_w, box_h, BLACK);

    DrawText("Default controls:", 20, 20, 10, BLACK);
    text_position += text_offset;

    DrawText("- Mouse Wheel to Zoom in-out", 40, text_position, 10, DARKGRAY);
    text_position += text_offset;

    DrawText("- Mouse Wheel Pressed to Pan", 40, text_position, 10, DARKGRAY);
    text_position += text_offset;

    DrawText("- Alt + Mouse Wheel Pressed to Rotate", 40, text_position, 10, DARKGRAY);
    text_position += text_offset;

    DrawText("- Alt + Ctrl + Mouse Wheel Pressed for Smooth Zoom", 40, text_position, 10, DARKGRAY);
    text_position += text_offset;

    DrawText("- Z to reset zoom", 40, text_position, 10, DARKGRAY);
    text_position += text_offset;

    DrawText("- G to only draw the grid lines", 40, text_position, 10, DARKGRAY);
    text_position += text_offset;

    DrawText("- L to enable or disable the grid lines", 40, text_position, 10, DARKGRAY);
    text_position += text_offset;

    DrawText("- Double click on a volume to show the AP", 40, text_position, 10, DARKGRAY);
    text_position += text_offset;

    DrawText("- R to restart simulation", 40, text_position, 12, BLACK);
    text_position += text_offset;

    DrawText("- Space to start or pause simulation", 40, text_position, 12, BLACK);
    text_position += text_offset;

    sprintf(tmp, "%lf ms", draw_config.time);
    if(draw_config.paused) {
        DrawText("Simulation paused:", 20, text_position, 16, BLACK);
    }
    else {
        DrawText("Simulation running:", 20, text_position, 16, BLACK);
    }
    DrawText(tmp, 170, text_position, 16, BLACK);
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

void init_and_open_visualization_window() {

    omp_set_lock(&draw_config.sleep_lock);

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);

    InitWindow(screenWidth, screenHeight, "Simulation visualization");

    font = LoadFont("misc/Roboto-Black.ttf");

    bool mesh_loaded = false;

    draw_config.grid_only = false;
    draw_config.grid_lines = true;

    selected_aps = NULL;
    hmdefault(selected_aps, NULL);

    Camera3D camera;
    camera.position = (Vector3){ 0.064882f, 0.165282f, 15.977825f };  // Camera position
    camera.target = (Vector3){ 0.137565f, 0.199405f, 0.181663f };
    camera.up       = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy     = 45.0f;                                // Camera field-of-view Y
    camera.type     = CAMERA_PERSPECTIVE;                  // Camera mode type

    SetCameraMode(camera, CAMERA_FREE); // Set a free camera mode

    SetTargetFPS(120);

    real_cpu scale = 1.0f;

    Ray ray;
    ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};

    Vector3 mesh_offset = (Vector3){ 0, 0, 0 };

    double mouse_timer = -1;

    while (!WindowShouldClose()) {

        if (IsKeyDown('Z')) {
            camera.target   = (Vector3){ 0.137565f, 0.199405f, 0.181663f };
        }

        UpdateCamera(&camera);

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
            mesh_loaded = false;
            ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
            ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
        }

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

        if(!draw_config.restart && draw_config.grid_to_draw && omp_test_lock(&draw_config.draw_lock)) {

            mesh_loaded = true;
            ClearBackground(GRAY);

            BeginMode3D(camera);

            if(!calc_center) {
                mesh_offset = find_mesh_center();
                scale = fmaxf(draw_config.grid_to_draw->side_length_x, fmaxf(draw_config.grid_to_draw->side_length_y, draw_config.grid_to_draw->side_length_z))/5.0f;
            }

            draw_alg_mesh(mesh_offset, scale, ray);
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
        else if(!mesh_loaded ){
            int posx = GetScreenWidth()/2-150;
            int posy = GetScreenHeight()/2-50;
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
