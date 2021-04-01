//
// Created by sachetto on 31/03/2021.
//

#ifndef MONOALG3D_C_RAYLIB_EXT_H
#define MONOALG3D_C_RAYLIB_EXT_H

#include "../3dparty/raylib/src/raymath.h"
#include <stdint.h>

void DrawCubeWithVisibilityMask(Vector3 position, float width, float height, float length, Color color, uint8_t visibility_mask);
void DrawCubeWiresWithVisibilityMask(Vector3 position, float width, float height, float length, Color color, uint8_t visibility_mask);
void DrawTextEx2(Font font, const char *text, Vector2 position, float fontSize, float spacing, Color tint);

#endif // MONOALG3D_C_RAYLIB_EXT_H
