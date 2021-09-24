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

void DrawMeshInstancedWithColors(Mesh mesh, Shader shader, Color* colors, Matrix *transforms, int grid_mask, int instances);

void DrawCubeWithVisibilityMask(Vector3 position, float width, float height, float length, Color color, uint8_t visibility_mask);
void DrawCubeWiresWithVisibilityMask(Vector3 position, float width, float height, float length, Color color, uint8_t visibility_mask);
void DrawTextEx2(Font font, const char *text, Vector2 position, float fontSize, float spacing, Color tint);

#endif // MONOALG3D_C_RAYLIB_EXT_H
