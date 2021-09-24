#include "raylib_ext.h"
#include "../3dparty/raylib/src/rlgl.h"
#include "../alg/cell/cell.h"

// Draw multiple mesh instances with material and different transforms
void DrawMeshInstancedWithColors(Mesh mesh, Shader shader, Color *colors, Matrix *transforms, int grid_mask, int instances)
{
#if defined(GRAPHICS_API_OPENGL_33) || defined(GRAPHICS_API_OPENGL_ES2)
    // Instancing required variables
    float16 *instanceTransforms = NULL;
    unsigned int instancesVboId = 0;

    float4 *colorsTransforms = NULL;
    unsigned int colorsVboId = 0;

    // Bind shader program
    rlEnableShader(shader.id);

    // Get a copy of current matrices to work with,
    // just in case stereo render is required and we need to modify them
    // NOTE: At this point the modelview matrix just contains the view matrix (camera)
    // That's because BeginMode3D() sets it and there is no model-drawing function
    // that modifies it, all use rlPushMatrix() and rlPopMatrix()
    Matrix matModel = MatrixIdentity();
    Matrix matView = rlGetMatrixModelview();
    Matrix matModelView = MatrixIdentity();
    Matrix matProjection = rlGetMatrixProjection();

    // Upload view and projection matrices (if locations available)
    if (shader.locs[SHADER_LOC_MATRIX_VIEW] != -1) rlSetUniformMatrix(shader.locs[SHADER_LOC_MATRIX_VIEW], matView);
    if (shader.locs[SHADER_LOC_MATRIX_PROJECTION] != -1) rlSetUniformMatrix(shader.locs[SHADER_LOC_MATRIX_PROJECTION], matProjection);

    // Create instances buffer
    instanceTransforms = (float16 *)RL_MALLOC(instances*sizeof(float16));
    colorsTransforms = (float4 *)RL_MALLOC(instances*sizeof(float4));

    // Fill buffer with instances transformations as float16 arrays
    for (int i = 0; i < instances; i++) instanceTransforms[i] = MatrixToFloatV(transforms[i]);
    for (int i = 0; i < instances; i++) {
        colorsTransforms[i].v[0] = colors[i].r/255.0;
        colorsTransforms[i].v[1] = colors[i].g/255.0;
        colorsTransforms[i].v[2] = colors[i].b/255.0;
        colorsTransforms[i].v[3] = colors[i].a;

    }

    // Enable mesh VAO to attach new buffer
    rlEnableVertexArray(mesh.vaoId);

    // This could alternatively use a static VBO and either glMapBuffer() or glBufferSubData().
    // It isn't clear which would be reliably faster in all cases and on all platforms,
    // anecdotally glMapBuffer() seems very slow (syncs) while glBufferSubData() seems
    // no faster, since we're transferring all the transform matrices anyway
    instancesVboId = rlLoadVertexBuffer(instanceTransforms, instances*sizeof(float16), false);

    // Instances transformation matrices are send to shader attribute location: SHADER_LOC_MATRIX_MODEL
    for (unsigned int i = 0; i < 4; i++)
    {
        rlEnableVertexAttribute(shader.locs[SHADER_LOC_MATRIX_MODEL] + i);
        rlSetVertexAttribute(shader.locs[SHADER_LOC_MATRIX_MODEL] + i, 4, RL_FLOAT, 0, sizeof(Matrix), (void *)(i*sizeof(Vector4)));
        rlSetVertexAttributeDivisor(shader.locs[SHADER_LOC_MATRIX_MODEL] + i, 1);
    }

    rlDisableVertexBuffer();
    rlDisableVertexArray();

    // Enable mesh VAO to attach new buffer
    rlEnableVertexArray(mesh.vaoId);
    colorsVboId = rlLoadVertexBuffer(colorsTransforms, instances*sizeof(float4), true);

    //Colors are send to shader attribute location: SHADER_LOC_VERTEX_COLOR
    rlEnableVertexAttribute(shader.locs[SHADER_LOC_VERTEX_COLOR]);
    rlSetVertexAttribute(shader.locs[SHADER_LOC_VERTEX_COLOR], 4, RL_FLOAT, 0, sizeof(float4), 0);
    rlSetVertexAttributeDivisor(shader.locs[SHADER_LOC_VERTEX_COLOR], 1);

    rlDisableVertexBuffer();
    rlDisableVertexArray();

    // Accumulate internal matrix transform (push/pop) and view matrix
    // NOTE: In this case, model instance transformation must be computed in the shader
    matModelView = MatrixMultiply(rlGetMatrixTransform(), matView);

    // Upload model normal matrix (if locations available)
    if (shader.locs[SHADER_LOC_MATRIX_NORMAL] != -1) rlSetUniformMatrix(shader.locs[SHADER_LOC_MATRIX_NORMAL], MatrixTranspose(MatrixInvert(matModel)));

    int dgrid_loc = GetShaderLocation(shader, "dgrid");
    if(dgrid_loc != 1) {
        rlSetUniform(dgrid_loc, (void *)&grid_mask, RL_SHADER_UNIFORM_INT, 1);
    }

    rlEnableVertexArray(mesh.vaoId);
    if (mesh.indices != NULL) rlEnableVertexBufferElement(mesh.vboId[6]);

    int eyeCount = 1;
    if (rlIsStereoRenderEnabled()) eyeCount = 2;

    for (int eye = 0; eye < eyeCount; eye++)
    {
        // Calculate model-view-projection matrix (MVP)
        Matrix matModelViewProjection = MatrixIdentity();
        if (eyeCount == 1) matModelViewProjection = MatrixMultiply(matModelView, matProjection);
        else
        {
            // Setup current eye viewport (half screen width)
            rlViewport(eye*rlGetFramebufferWidth()/2, 0, rlGetFramebufferWidth()/2, rlGetFramebufferHeight());
            matModelViewProjection = MatrixMultiply(MatrixMultiply(matModelView, rlGetMatrixViewOffsetStereo(eye)), rlGetMatrixProjectionStereo(eye));
        }

        // Send combined model-view-projection matrix to shader
        rlSetUniformMatrix(shader.locs[SHADER_LOC_MATRIX_MVP], matModelViewProjection);

        // Draw mesh instanced
        if (mesh.indices != NULL) {
            rlDrawVertexArrayElementsInstanced(0, mesh.triangleCount*3, 0, instances);
        }
        else {
            rlDrawVertexArrayInstanced(0, mesh.vertexCount, instances);
        }
    }

    // Disable all possible vertex array objects (or VBOs)
    rlDisableVertexArray();
    rlDisableVertexBuffer();
    rlDisableVertexBufferElement();

    // Disable shader program
    rlDisableShader();

    // Remove instance transforms buffer
    rlUnloadVertexBuffer(instancesVboId);
    rlUnloadVertexBuffer(colorsVboId);
    RL_FREE(instanceTransforms);
    RL_FREE(colorsTransforms);
#endif
}



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
/*
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
int codepoint = GetCodepoint(&text[i], &codepointByteCount);
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
Rectangle rec = { position.x + textOffsetX + font.recs[index].offsetX*scaleFactor,
position.y - textOffsetY - font.recs[index].offsetY*scaleFactor,
font.recs[index].width*scaleFactor,
font.recs[index].height*scaleFactor };

DrawTexturePro(font.texture, font.recs[index], rec, (Vector2){ 0, 0 }, -90.0f, tint);
}

if (font.rect[index].advanceX == 0) textOffsetY += (int)(font.recs[index].width*scaleFactor + spacing);
else textOffsetY += (int)((float)font.chars[index].advanceX*scaleFactor + spacing);
}

i += (codepointByteCount - 1);   // Move text bytes counter to next codepoint
}
}
*/

