//
// Created by sachetto on 31/03/2021.
//

#ifndef MONOALG3D_C_RAYLIB_EXT_H
#define MONOALG3D_C_RAYLIB_EXT_H

#include "../3dparty/raylib/src/raylib.h"
#include "../3dparty/raylib/src/raymath.h"
#include <stdint.h>

struct draw_context {
    Shader shader;
    Mesh mesh;
    Color *colors;
    float16 *instance_transforms;
    unsigned int instances_vbo;
    unsigned int colors_vbo;
    int instance_capacity;
    int grid_mask_location;
};

void UnloadMeshInstanceBuffers(struct draw_context *draw_context);

void DrawMeshInstancedWithColors(struct draw_context *draw_context, int grid_mask, int instances);

#endif // MONOALG3D_C_RAYLIB_EXT_H
