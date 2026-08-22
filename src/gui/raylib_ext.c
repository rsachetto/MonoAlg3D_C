#include "raylib_ext.h"
#include "../3dparty/raylib/src/rlgl.h"
#include "../alg/cell/cell.h"

// Draw multiple mesh instances with different transforms and colors
void DrawMeshInstancedWithColors(struct draw_context *draw_context, int grid_mask, int instances) {

    if(instances <= 0) return;

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

    // Colors may change every frame. Geometry generally does not, so avoid the
    // matrix conversion and the much larger GPU upload until the mesh changes.
    for(int i = 0; i < instances; i++) {
        if(draw_context->transforms_dirty) {
            draw_context->instance_transforms[i] = MatrixToFloatV(draw_context->translations[i]);
        }
        draw_context->colors_transforms[i].v[0] = (float)draw_context->colors[i].r / 255.0f;
        draw_context->colors_transforms[i].v[1] = (float)draw_context->colors[i].g / 255.0f;
        draw_context->colors_transforms[i].v[2] = (float)draw_context->colors[i].b / 255.0f;
        draw_context->colors_transforms[i].v[3] = (float)draw_context->colors[i].a / 255.0f;
    }

    // Enable mesh VAO to attach new buffer
    rlEnableVertexArray(draw_context->mesh.vaoId);

    // Keep the instance buffers alive between frames. Recreating them here used to
    // force two GPU allocations and two deletions for every rendered frame, which is
    // particularly expensive for large meshes.
    if(draw_context->instances_vbo_id == 0 || instances > draw_context->vbo_capacity) {
        UnloadDrawContextBuffers(draw_context);
        draw_context->instances_vbo_id =
            rlLoadVertexBuffer(draw_context->instance_transforms, (int)(instances * sizeof(float16)), true);
        draw_context->colors_vbo_id =
            rlLoadVertexBuffer(draw_context->colors_transforms, (int)(instances * sizeof(float4)), true);
        draw_context->vbo_capacity = instances;
    } else {
        if(draw_context->transforms_dirty) {
            rlUpdateVertexBuffer(draw_context->instances_vbo_id, draw_context->instance_transforms,
                                 (int)(instances * sizeof(float16)), 0);
        }
        rlUpdateVertexBuffer(draw_context->colors_vbo_id, draw_context->colors_transforms,
                             (int)(instances * sizeof(float4)), 0);
    }

    rlEnableVertexBuffer(draw_context->instances_vbo_id);

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
    rlEnableVertexBuffer(draw_context->colors_vbo_id);

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

    draw_context->transforms_dirty = false;

}

void UnloadDrawContextBuffers(struct draw_context *draw_context) {
    if(draw_context->instances_vbo_id != 0) {
        rlUnloadVertexBuffer(draw_context->instances_vbo_id);
        draw_context->instances_vbo_id = 0;
    }

    if(draw_context->colors_vbo_id != 0) {
        rlUnloadVertexBuffer(draw_context->colors_vbo_id);
        draw_context->colors_vbo_id = 0;
    }

    draw_context->vbo_capacity = 0;
}
