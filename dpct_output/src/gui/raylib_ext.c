#include "raylib_ext.h"
#include "../3dparty/raylib/src/rlgl.h"
#include "../alg/cell/cell.h"

// Draw multiple mesh instances with different transforms and colors
void DrawMeshInstancedWithColors(struct draw_context *draw_context, int grid_mask, int instances) {

    // Instancing required variables
    unsigned int instancesVboId;
    unsigned int colorsVboId;

    // Bind shader program
    rlEnableShader(draw_context->shader.id);

    // Get a copy of current matrices to work with,
    // just in case stereo render is required and we need to modify them
    // NOTE: At this point the modelview matrix just contains the view matrix (camera)
    // That's because BeginMode3D() sets it and there is no model-drawing function
    // that modifies it, all use rlPushMatrix() and rlPopMatrix()
    static const Matrix matModel = {1.0f, 0.0f, 0.0f, 0.0f,
                                    0.0f, 1.0f, 0.0f, 0.0f,
                                    0.0f, 0.0f, 1.0f, 0.0f,
                                    0.0f, 0.0f, 0.0f, 1.0f};

    Matrix matView = rlGetMatrixModelview();
    Matrix matModelView;
    Matrix matProjection = rlGetMatrixProjection();

    rlSetUniformMatrix(draw_context->shader.locs[SHADER_LOC_MATRIX_VIEW], matView);
    rlSetUniformMatrix(draw_context->shader.locs[SHADER_LOC_MATRIX_PROJECTION], matProjection);

    // Fill buffer with instances transformations as float16 arrays
    for(int i = 0; i < instances; i++) {
        draw_context->instance_transforms[i] = MatrixToFloatV(draw_context->translations[i]);
        draw_context->colors_transforms[i].v[0] = (float)draw_context->colors[i].r / 255.0f;
        draw_context->colors_transforms[i].v[1] = (float)draw_context->colors[i].g / 255.0f;
        draw_context->colors_transforms[i].v[2] = (float)draw_context->colors[i].b / 255.0f;
        draw_context->colors_transforms[i].v[3] = (float)draw_context->colors[i].a / 255.0f;
    }

    // Enable mesh VAO to attach new buffer
    rlEnableVertexArray(draw_context->mesh.vaoId);

    // This could alternatively use a static VBO and either glMapBuffer() or glBufferSubData().
    // It isn't clear which would be reliably faster in all cases and on all platforms,
    // anecdotally glMapBuffer() seems very slow (syncs) while glBufferSubData() seems
    // no faster, since we're transferring all the transform matrices anyway
    instancesVboId = rlLoadVertexBuffer(draw_context->instance_transforms, (int)(instances * sizeof(float16)), false);

    // Instances transformation matrices are send to shader attribute location: SHADER_LOC_MATRIX_MODEL
    for(unsigned int i = 0; i < 4; i++) {
        rlEnableVertexAttribute(draw_context->shader.locs[SHADER_LOC_MATRIX_MODEL] + i);
        rlSetVertexAttribute(draw_context->shader.locs[SHADER_LOC_MATRIX_MODEL] + i, 4, RL_FLOAT, 0, sizeof(Matrix), (void *)(i * sizeof(Vector4)));
        rlSetVertexAttributeDivisor(draw_context->shader.locs[SHADER_LOC_MATRIX_MODEL] + i, 1);
    }

    rlDisableVertexBuffer();
    rlDisableVertexArray();

    // Enable mesh VAO to attach new buffer
    rlEnableVertexArray(draw_context->mesh.vaoId);
    colorsVboId = rlLoadVertexBuffer(draw_context->colors_transforms, (int)(instances * sizeof(float4)), true);

    // Colors are send to shader attribute location: SHADER_LOC_VERTEX_COLOR
    rlEnableVertexAttribute(draw_context->shader.locs[SHADER_LOC_VERTEX_COLOR]);
    rlSetVertexAttribute(draw_context->shader.locs[SHADER_LOC_VERTEX_COLOR], 4, RL_FLOAT, 0, sizeof(float4), 0);
    rlSetVertexAttributeDivisor(draw_context->shader.locs[SHADER_LOC_VERTEX_COLOR], 1);

    rlDisableVertexBuffer();
    rlDisableVertexArray();

    // Accumulate internal matrix transform (push/pop) and view matrix
    // NOTE: In this case, model instance transformation must be computed in the shader
    matModelView = MatrixMultiply(rlGetMatrixTransform(), matView);
    rlSetUniformMatrix(draw_context->shader.locs[SHADER_LOC_MATRIX_NORMAL], matModel);

    int dgrid_loc = GetShaderLocation(draw_context->shader, "dgrid");
    rlSetUniform(dgrid_loc, (void *)&grid_mask, RL_SHADER_UNIFORM_INT, 1);

    rlEnableVertexArray(draw_context->mesh.vaoId);
    rlEnableVertexBufferElement(draw_context->mesh.vboId[6]);

    // Calculate model-view-projection matrix (MVP)
    Matrix matModelViewProjection;
    matModelViewProjection = MatrixMultiply(matModelView, matProjection);

    // Send combined model-view-projection matrix to shader
    rlSetUniformMatrix(draw_context->shader.locs[SHADER_LOC_MATRIX_MVP], matModelViewProjection);
    rlDrawVertexArrayElementsInstanced(0, draw_context->mesh.triangleCount * 3, 0, instances);

    // Disable all possible vertex array objects (or VBOs)
    rlDisableVertexArray();
    rlDisableVertexBuffer();
    rlDisableVertexBufferElement();

    // Disable shader program
    rlDisableShader();

    // Remove instance transforms buffer
    rlUnloadVertexBuffer(instancesVboId);
    rlUnloadVertexBuffer(colorsVboId);

}
