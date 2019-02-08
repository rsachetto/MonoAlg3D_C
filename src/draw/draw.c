//
// Created by sachetto on 11/11/17.
//

#include <float.h>
#include <GL/gl.h>
#include "draw.h"
#include "../raylib/raylib.h"
#include "../raylib/rlgl.h"

bool calc_center = false;

Color get_color(double value)
{
//    #define NUM_COLORS 4
//    static float color[NUM_COLORS][3] = { {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };

    int idx1;        // |-- Our desired color will be between these two indexes in "color".
    int idx2;        // |
    double fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

    if(value <= 0)      {  idx1 = idx2 = 0;            }    // accounts for an input <=0
    else if(value >= 1)  {  idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
    else
    {
        value = value * (NUM_COLORS-1);        // Will multiply value by 3.
        idx1  = (int)floor(value);                  // Our desired color will be after this index.
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

static Vector3 find_mesh_center() {

    struct grid *grid_to_draw = draw_config.grid_to_draw;

    uint32_t n_active = grid_to_draw->num_active_cells;
    struct cell_node **ac = grid_to_draw->active_cells;
    struct cell_node *grid_cell;

    double max_x, max_y, max_z;
    double min_x, min_y, min_z;

    max_x = FLT_MIN;
    max_y = FLT_MIN;
    max_z = FLT_MIN;

    min_x = FLT_MAX;
    min_y = FLT_MAX;
    min_z = FLT_MAX;

    Vector3 result = (Vector3){0.0, 0.0, 0.0};

    if (grid_to_draw) {

        calc_center = true;

        if (ac) {
//            #pragma omp parallel for
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

        printf("%lf, %lf, %lf", result.x, result.y, result.z);
        calc_center = true;

    }



    return result;

}

void draw_coordinates(){

    // draw some lines
    rlColor3f(1.0,0.0,0.0); // red x
    rlBegin(GL_LINES);
    // x aix

    rlVertex3f(-4.0f, 0.0f, 0.0f);
    rlVertex3f(4.0f, 0.0f, 0.0f);

    rlVertex3f(4.0f, 0.0f, 0.0f);
    rlVertex3f(3.0f, 1.0f, 0.0f);

    rlVertex3f(4.0f, 0.0f, 0.0f);
    rlVertex3f(3.0f, -1.0f, 0.0f);
    rlEnd();

    // y 
    rlColor3f(0.0f,1.0f,0.0f); // green y
    rlBegin(GL_LINES);
    rlVertex3f(0.0f, -4.0f, 0.0f);
    rlVertex3f(0.0f, 4.0f, 0.0f);

    rlVertex3f(0.0, 4.0f, 0.0f);
    rlVertex3f(1.0, 3.0f, 0.0f);

    rlVertex3f(0.0, 4.0f, 0.0f);
    rlVertex3f(-1.0, 3.0f, 0.0f);
    rlEnd();

    // z 
    rlColor3f(0.0,0.0,1.0); // blue z
    rlBegin(GL_LINES);
    rlVertex3f(0.0, 0.0f ,-4.0f );
    rlVertex3f(0.0, 0.0f ,4.0f );


    rlVertex3f(0.0, 0.0f ,4.0f );
    rlVertex3f(0.0, 1.0f ,3.0f );

    rlVertex3f(0.0, 0.0f ,4.0f );
    rlVertex3f(0.0, -1.0f ,3.0f );
    rlEnd();

}


static void draw_alg_mesh(Vector3 mesh_offset) {

    struct grid *grid_to_draw = draw_config.grid_to_draw;

   // draw_coordinates();

    //return;

    if (grid_to_draw) {

        uint32_t n_active = grid_to_draw->num_active_cells;
        struct cell_node **ac = grid_to_draw->active_cells;
        struct cell_node *grid_cell;

//        double last_non_nan = 0.0;

        if (ac) {
            //#pragma omp parallel for
            for (int i = 0; i < n_active; i++) {
                grid_cell = ac[i];

                Vector3 cubePosition;

                cubePosition.x = (grid_cell->center_x - mesh_offset.x)/2000.0f;
                cubePosition.y = (grid_cell->center_y - mesh_offset.y)/2000.0f;
                cubePosition.z = (grid_cell->center_z - mesh_offset.z)/2000.0f;

                Vector3 cubeSize;
                cubeSize.x = grid_cell->dx/2000.0f;
                cubeSize.y = grid_cell->dy/2000.0f;
                cubeSize.z = grid_cell->dz/2000.0f;

                double v = (grid_cell->v - draw_config.min_v)/(draw_config.max_v - draw_config.min_v);

//                if(!isnan(v)) {
//                    last_non_nan = v;
//                }
//                else {
//                    v = last_non_nan;
//                }

                Color color = get_color(v);

                if(draw_config.grid_only) {
                    DrawCubeWiresV(cubePosition, cubeSize, color);
                }
                else {
                    DrawCubeV(cubePosition, cubeSize, color);
                    if(draw_config.grid_lines) {
                        DrawCubeWiresV(cubePosition, cubeSize, BLACK);
                    }
                }


            }


        }
    }

}

void init_opengl() {
    // Initialization
    //--------------------------------------------------------------------------------------

    int screenWidth = 1280;
    int screenHeight = 720;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);

    InitWindow(screenWidth, screenHeight, "Simulation visualization");

    draw_config.grid_only = false;
    draw_config.grid_lines = true;

// Define the camera to look into our 3d world
    Camera3D camera;
    camera.position = (Vector3){ 0.064882f, 0.165282f, 15.977825f };  // Camera position
    camera.target   = (Vector3){ 0.137565f, 0.199405f, 0.181663f };      // Camera looking at point
    camera.up       = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy     = 45.0f;                                // Camera field-of-view Y
    camera.type     = CAMERA_PERSPECTIVE;                  // Camera mode type


    SetCameraMode(camera, CAMERA_FREE); // Set a free camera mode

    SetTargetFPS(60);   // Set our game to run at 60 frames-per-second

    Vector3 mesh_offset = (Vector3){ 0.0f, 0.0f, 0.0f };
    //--------------------------------------------------------------------------------------
    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        UpdateCamera(&camera);          // Update camera

        if (IsKeyDown('Z')) camera.target = (Vector3){ 3.5f, -2.0f, -1.0f };

        if (IsKeyPressed('G')) draw_config.grid_only = !draw_config.grid_only;

        if (IsKeyPressed('L')) draw_config.grid_lines = !draw_config.grid_lines;

        // Draw
        //----------------------------------------------------------------------------------

//        printf("POSITION: %lf, %lf, %lf\n", camera.position.x, camera.position.y, camera.position.z);
//        printf("TARGET: %lf, %lf, %lf\n", camera.target.x, camera.target.y, camera.target.z);

        BeginDrawing();

        if(draw_config.grid_to_draw && omp_test_lock(&draw_config.draw_lock)) {

            BeginMode3D(camera);
           // DrawGrid(10, 1.0f);
            ClearBackground(GRAY);
            if(!calc_center) {
                mesh_offset = find_mesh_center();
            }
            draw_alg_mesh(mesh_offset);
            omp_unset_lock(&draw_config.draw_lock);


            EndMode3D();

            //DrawFPS(10, 360);

            DrawRectangle(10, 10, 320, 173, Fade(SKYBLUE, 0.5f));
            DrawRectangleLines(10, 10, 320, 173, BLUE);

            DrawText("Free camera default controls:", 20, 20, 10, BLACK);
            DrawText("- Mouse Wheel to Zoom in-out", 40, 40, 10, DARKGRAY);
            DrawText("- Mouse Wheel Pressed to Pan", 40, 60, 10, DARKGRAY);
            DrawText("- Alt + Mouse Wheel Pressed to Rotate", 40, 80, 10, DARKGRAY);
            DrawText("- Alt + Ctrl + Mouse Wheel Pressed for Smooth Zoom", 40, 100, 10, DARKGRAY);
            DrawText("- Z to reset zoom", 40, 120, 10, DARKGRAY);
            DrawText("- G to only draw the grid lines", 40, 140, 10, DARKGRAY);
            DrawText("- L to enable or disable the grid lines", 40, 160, 10, DARKGRAY);
        }

        EndDrawing();
    }
        //----------------------------------------------------------------------------------


    // De-Initialization
    omp_destroy_lock(&draw_config.draw_lock);


    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------



}