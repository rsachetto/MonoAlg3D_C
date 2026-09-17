#include "raylib_ext.h"
#include "../3dparty/raylib/src/rlgl.h"
#include <limits.h>
#include <stddef.h>

void UnloadMeshInstanceBuffers(struct draw_context *draw_context) {
    if(draw_context->instances_vbo) rlUnloadVertexBuffer(draw_context->instances_vbo);
    if(draw_context->colors_vbo) rlUnloadVertexBuffer(draw_context->colors_vbo);
    draw_context->instances_vbo = 0;
    draw_context->colors_vbo = 0;
    draw_context->instance_capacity = 0;
}

// Draw multiple mesh instances with different transforms and colors
void DrawMeshInstancedWithColors(struct draw_context *draw_context, int grid_mask, int instances) {

    if(instances <= 0) return;
    // rlgl uses signed int byte counts for buffer allocation and updates.
    if(instances > INT_MAX / (int)sizeof(float16)) return;

    if(instances > draw_context->instance_capacity) {
        UnloadMeshInstanceBuffers(draw_context);
        draw_context->instances_vbo = rlLoadVertexBuffer(NULL, instances * sizeof(float16), true);
        draw_context->colors_vbo = rlLoadVertexBuffer(NULL, instances * sizeof(Color), true);
        if(!draw_context->instances_vbo || !draw_context->colors_vbo) {
            UnloadMeshInstanceBuffers(draw_context);
            rlDisableVertexBuffer();
            return;
        }
        draw_context->instance_capacity = instances;
    }

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

    // Enable mesh VAO to attach new buffer
    rlEnableVertexArray(draw_context->mesh.vaoId);

    rlUpdateVertexBuffer(draw_context->instances_vbo, draw_context->instance_transforms,
                         instances * sizeof(float16), 0);

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
    rlUpdateVertexBuffer(draw_context->colors_vbo, draw_context->colors, instances * sizeof(Color), 0);

    // Colors are send to shader attribute location: SHADER_LOC_VERTEX_COLOR
    rlEnableVertexAttribute(draw_context->shader.locs[SHADER_LOC_VERTEX_COLOR]);
    rlSetVertexAttribute(draw_context->shader.locs[SHADER_LOC_VERTEX_COLOR], 4, RL_UNSIGNED_BYTE, true, sizeof(Color), 0);
    rlSetVertexAttributeDivisor(draw_context->shader.locs[SHADER_LOC_VERTEX_COLOR], 1);

    rlDisableVertexBuffer();
    rlDisableVertexArray();

    // Accumulate internal matrix transform (push/pop) and view matrix
    // NOTE: In this case, model instance transformation must be computed in the shader
    matModelView = MatrixMultiply(rlGetMatrixTransform(), matView);
    rlSetUniformMatrix(draw_context->shader.locs[SHADER_LOC_MATRIX_NORMAL], matModel);

    rlSetUniform(draw_context->grid_mask_location, (void *)&grid_mask, RL_SHADER_UNIFORM_INT, 1);

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


}
