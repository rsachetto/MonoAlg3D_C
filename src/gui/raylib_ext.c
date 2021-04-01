#include "raylib_ext.h"
#include "../3dparty/raylib/src/rlgl.h"
#include "../alg/cell/cell.h"

void DrawCubeWithVisibilityMask(Vector3 position, float width, float height, float length, Color color, uint8_t visibility_mask) {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;

    rlCheckRenderBatchLimit(36);

    rlPushMatrix();
    // NOTE: Transformation is applied in inverse order (scale -> rotate -> translate)
    rlTranslatef(position.x, position.y, position.z);

    rlBegin(RL_TRIANGLES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    if(visibility_mask & FRONT_IS_VISIBLE) {
        // Front face
        rlVertex3f(x - width/2, y - height/2, z + length/2);  // Bottom Left
        rlVertex3f(x + width/2, y - height/2, z + length/2);  // Bottom Right
        rlVertex3f(x - width/2, y + height/2, z + length/2);  // Top Left

        rlVertex3f(x + width/2, y + height/2, z + length/2);  // Top Right
        rlVertex3f(x - width/2, y + height/2, z + length/2);  // Top Left
        rlVertex3f(x + width/2, y - height/2, z + length/2);  // Bottom Right
    }

    if(visibility_mask & BACK_IS_VISIBLE) {
        // Back face
        rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Bottom Left
        rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left
        rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Bottom Right

        rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right
        rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Bottom Right
        rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left
    }

    if(visibility_mask & TOP_IS_VISIBLE) {
        // Top face
        rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left
        rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Bottom Left
        rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Bottom Right

        rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right
        rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left
        rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Bottom Right
    }

    if(visibility_mask & DOWN_IS_VISIBLE) {
        // Bottom face
        rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Top Left
        rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Bottom Right
        rlVertex3f(x - width / 2, y - height / 2, z + length / 2); // Bottom Left

        rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Top Right
        rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Bottom Right
        rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Top Left
    }

    if(visibility_mask & RIGHT_IS_VISIBLE) {
        // Right face
        rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Bottom Right
        rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right
        rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Top Left

        rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Bottom Left
        rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Bottom Right
        rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Top Left
    }

    if(visibility_mask & LEFT_IS_VISIBLE) {
        // Left face
        rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Bottom Right
        rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Top Left
        rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Right

        rlVertex3f(x - width / 2, y - height / 2, z + length / 2); // Bottom Left
        rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Top Left
        rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Bottom Right
    }
    rlEnd();
    rlPopMatrix();
}

void DrawCubeWiresWithVisibilityMask(Vector3 position, float width, float height, float length, Color color, uint8_t visibility_mask) {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;

    rlCheckRenderBatchLimit(36);

    rlPushMatrix();
    rlTranslatef(position.x, position.y, position.z);

    rlBegin(RL_LINES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    bool front_visible = visibility_mask & FRONT_IS_VISIBLE;
    bool back_visible  = visibility_mask & BACK_IS_VISIBLE;
    bool top_visible   = visibility_mask & TOP_IS_VISIBLE;
    bool down_visible  = visibility_mask & DOWN_IS_VISIBLE;
    bool left_visible  = visibility_mask & LEFT_IS_VISIBLE;
    bool right_visible = visibility_mask & RIGHT_IS_VISIBLE;

    if(front_visible) {
        // Front Face -----------------------------------------------------
        // Bottom Line
        rlVertex3f(x - width / 2, y - height / 2, z + length / 2); // Bottom Left
        rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Bottom Right

        // Left Line
        rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Bottom Right
        rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Top Right

        // Top Line
        rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Top Right
        rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Top Left

        // Right Line
        rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Top Left
        rlVertex3f(x - width / 2, y - height / 2, z + length / 2); // Bottom Left
    }

    if(back_visible) {
        // Back Face ------------------------------------------------------
        // Bottom Line
        rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Bottom Left
        rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Bottom Right

        // Left Line
        rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Bottom Right
        rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right

        // Top Line
        rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right
        rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left

        // Right Line
        rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left
        rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Bottom Left
    }

    if(top_visible) {

        if(!front_visible) {
            // Top Line
            rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Top Right
            rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Top Left
        }

        if(!back_visible) {
            // Top Line
            rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right
            rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left
        }

        // Top Face -------------------------------------------------------
        // Left Line
        rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Top Left Front
        rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left Back

        // Right Line
        rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Top Right Front
        rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right Back
    }
    if(down_visible) {

        if(!front_visible) {
            // Bottom Line
            rlVertex3f(x - width / 2, y - height / 2, z + length / 2); // Bottom Left
            rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Bottom Right
        }

        if(!back_visible) {
            // Bottom Line
            rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Bottom Left
            rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Bottom Right
        }

        // Bottom Face  ---------------------------------------------------
        // Left Line
        rlVertex3f(x - width / 2, y - height / 2, z + length / 2); // Top Left Front
        rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Top Left Back

        // Right Line
        rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Top Right Front
        rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Top Right Back
    }

    if(right_visible || left_visible) {
        if(!front_visible) {
            // Left Line
            rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Bottom Right
            rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Top Right

            // Right Line
            rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Top Left
            rlVertex3f(x - width / 2, y - height / 2, z + length / 2); // Bottom Left
        }

        if(!back_visible) {

            // Left Line
            rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Bottom Right
            rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right

            // Right Line
            rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left
            rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Bottom Left
        }

        if(!top_visible) {
            // Left Line
            rlVertex3f(x - width / 2, y + height / 2, z + length / 2); // Top Left Front
            rlVertex3f(x - width / 2, y + height / 2, z - length / 2); // Top Left Back

            // Right Line
            rlVertex3f(x + width / 2, y + height / 2, z + length / 2); // Top Right Front
            rlVertex3f(x + width / 2, y + height / 2, z - length / 2); // Top Right Back
        }
        if(!down_visible) {
            // Bottom Face  ---------------------------------------------------
            // Left Line
            rlVertex3f(x - width / 2, y - height / 2, z + length / 2); // Top Left Front
            rlVertex3f(x - width / 2, y - height / 2, z - length / 2); // Top Left Back

            // Right Line
            rlVertex3f(x + width / 2, y - height / 2, z + length / 2); // Top Right Front
            rlVertex3f(x + width / 2, y - height / 2, z - length / 2); // Top Right Back
        }

    }
    rlEnd();
    rlPopMatrix();
}

void DrawTextEx2(Font font, const char *text, Vector2 position, float fontSize, float spacing, Color tint)
{
    unsigned int length = TextLength(text);      // Total length in bytes of the text, scanned by codepoints in loop

    int textOffsetY = 0;            // Offset between lines (on line break '\n')
    float textOffsetX = 0.0f;       // Offset X to next character to draw

    float scaleFactor = fontSize/font.baseSize;     // Character quad scaling factor

    for (int i = 0; i < length; i++)
    {
        // Get next codepoint from byte sds and glyph index in font
        int codepointByteCount = 0;
        int codepoint = GetNextCodepoint(&text[i], &codepointByteCount);
        int index = GetGlyphIndex(font, codepoint);

        // NOTE: Normally we exit the decoding sequence as soon as a bad byte is found (and return 0x3f)
        // but we need to draw all of the bad bytes using the '?' symbol moving one byte
        if (codepoint == 0x3f) codepointByteCount = 1;

        if (codepoint == '\n')
        {
            textOffsetY += (int)((font.baseSize + font.baseSize/2)*scaleFactor);
            textOffsetX = 0.0f;
        }
        else
        {
            if ((codepoint != ' ') && (codepoint != '\t'))
            {
                Rectangle rec = { position.x + textOffsetX + font.chars[index].offsetX*scaleFactor,
                                  position.y - textOffsetY - font.chars[index].offsetY*scaleFactor,
                                  font.recs[index].width*scaleFactor,
                                  font.recs[index].height*scaleFactor };

                DrawTexturePro(font.texture, font.recs[index], rec, (Vector2){ 0, 0 }, -90.0f, tint);
            }

            if (font.chars[index].advanceX == 0) textOffsetY += (int)(font.recs[index].width*scaleFactor + spacing);
            else textOffsetY += (int)((float)font.chars[index].advanceX*scaleFactor + spacing);
        }

        i += (codepointByteCount - 1);   // Move text bytes counter to next codepoint
    }
}