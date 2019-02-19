//
// Created by sachetto on 11/11/17.
//

#include <float.h>
#include "draw.h"


#include "../raylib/src/raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "../raylib/src/raygui.h"

bool calc_center = false;
Vector3 *cube_positions = NULL;

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



static void draw_alg_mesh(Vector3 mesh_offset, float scale, Ray ray) {

    struct grid *grid_to_draw = draw_config.grid_to_draw;

    Vector3 cubePosition;
    Vector3 cubeSize;
    Color color;

    double max_v = draw_config.max_v;
    double min_v = draw_config.min_v;

    bool grid_only = draw_config.grid_only;
    bool grid_lines = draw_config.grid_lines;

    bool collision = false;

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

                cubePosition.x = (grid_cell->center_x - mesh_offset.x)/scale;
                cubePosition.y = (grid_cell->center_y - mesh_offset.y)/scale;
                cubePosition.z = (grid_cell->center_z - mesh_offset.z)/scale;

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
                    }
                }
            }
        }

    }
}

const int screenWidth = 1280;
const int screenHeight = 720;

void init_opengl() {

    SetConfigFlags(FLAG_WINDOW_RESIZABLE);

    InitWindow(screenWidth, screenHeight, "Simulation visualization");

    draw_config.grid_only = false;
    draw_config.grid_lines = true;

    Camera3D camera;
    camera.position = (Vector3){ 0.064882f, 0.165282f, 15.977825f };  // Camera position
    camera.target   = (Vector3){ 0.137565f, 0.199405f, 0.181663f };      // Camera looking at point
    camera.up       = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy     = 45.0f;                                // Camera field-of-view Y
    camera.type     = CAMERA_PERSPECTIVE;                  // Camera mode type

    SetCameraMode(camera, CAMERA_FREE); // Set a free camera mode

    SetTargetFPS(120);

    Vector3 mesh_offset = (Vector3){ 0.0f, 0.0f, 0.0f };
    float scale = 1.0f;

    char tmp[100];
    bool mesh_loaded = false;

    Ray ray = {FLT_MAX, FLT_MAX, FLT_MAX};        // Picking line ray

    while (!WindowShouldClose())    // Detect window close button or ESC key
    {

        if (IsKeyDown('Z')) {
            camera.target   = (Vector3){ 0.137565f, 0.199405f, 0.181663f };
        }

        UpdateCamera(&camera);

        if (IsKeyPressed('G')) draw_config.grid_only = !draw_config.grid_only;

        if (IsKeyPressed('L')) draw_config.grid_lines = !draw_config.grid_lines;

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            ray = GetMouseRay(GetMousePosition(), camera);
        }

        if(draw_config.grid_to_draw && omp_test_lock(&draw_config.draw_lock)) {

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

                DrawRectangle(10, 10, 320, 193, SKYBLUE);
                DrawRectangleLines(10, 10, 320, 193, BLUE);
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

            else {

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
        }
        else if(!mesh_loaded){
            int posx = GetScreenWidth()/2-150;
            int posy = GetScreenHeight()/2-50;
            ClearBackground(GRAY);
            DrawRectangle(posx, posy, 320, 20, SKYBLUE);
            DrawRectangleLines(posx, posy, 320, 20, BLUE);
            DrawText("Loading Mesh...", posx+80, posy, 20, BLACK);
        }

        if(draw_config.simulating) {
            if(!draw_config.paused) {
                if (GuiButton((Rectangle){ GetScreenWidth() - 170, GetScreenHeight() - 50, 150, 30 }, "Pause Simulation")) {
                    draw_config.paused = true;
                    printf("PAUSING\n");
                }
            } else {
                if (GuiButton((Rectangle){ GetScreenWidth() - 170, GetScreenHeight() - 50, 150, 30 }, "Continue Simulation")) {
                    draw_config.paused = false;
                }
            }
            GuiEnable();
        }

        EndDrawing();

    }

    omp_destroy_lock(&draw_config.draw_lock);

    CloseWindow();

}
