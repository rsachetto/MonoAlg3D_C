//
// Created by sachetto on 31/03/2021.
//

#ifndef MONOALG3D_C_RAYLIB_EXT_H
#define MONOALG3D_C_RAYLIB_EXT_H

#include "../3dparty/raylib/src/raylib.h"
#include "../3dparty/raylib/src/raymath.h"
#include <stdint.h>

typedef struct float4 {
    float v[4];
} float4;

struct draw_context {
    Shader shader;
    Mesh mesh;
    Matrix *translations;
    Color *colors;
    float16 *instance_transforms;
    float4 *colors_transforms;
};

void DrawMeshInstancedWithColors(struct draw_context *draw_context, int grid_mask, int instances);

#endif // MONOALG3D_C_RAYLIB_EXT_H
