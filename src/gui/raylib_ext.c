#include "raylib_ext.h"
#include "../3dparty/raylib/src/rlgl.h"
#include "../alg/cell/cell.h"

// Draw multiple mesh instances with different transforms and colors
void DrawMeshInstancedWithColors(Mesh mesh, Shader shader, Color *colors, Matrix *transforms, int grid_mask, int instances) {
#if defined(GRAPHICS_API_OPENGL_33) || defined(GRAPHICS_API_OPENGL_ES2)
    // Instancing required variables
    float16 *instanceTransforms = NULL;
    unsigned int instancesVboId;

    float4 *colorsTransforms = NULL;
    unsigned int colorsVboId;

    // Bind shader program
    rlEnableShader(shader.id);

    // Get a copy of current matrices to work with,
    // just in case stereo render is required and we need to modify them
    // NOTE: At this point the modelview matrix just contains the view matrix (camera)
    // That's because BeginMode3D() sets it and there is no model-drawing function
    // that modifies it, all use rlPushMatrix() and rlPopMatrix()
    Matrix matModel = MatrixIdentity();
    Matrix matView = rlGetMatrixModelview();
    Matrix matModelView;
    Matrix matProjection = rlGetMatrixProjection();

    // Upload view and projection matrices (if locations available)
    if(shader.locs[SHADER_LOC_MATRIX_VIEW] != -1)
        rlSetUniformMatrix(shader.locs[SHADER_LOC_MATRIX_VIEW], matView);
    if(shader.locs[SHADER_LOC_MATRIX_PROJECTION] != -1)
        rlSetUniformMatrix(shader.locs[SHADER_LOC_MATRIX_PROJECTION], matProjection);

    // Create instances buffer
    instanceTransforms = (float16 *)RL_MALLOC(instances * sizeof(float16));
    colorsTransforms = (float4 *)RL_MALLOC(instances * sizeof(float4));

    // Fill buffer with instances transformations as float16 arrays
    for(int i = 0; i < instances; i++)
        instanceTransforms[i] = MatrixToFloatV(transforms[i]);
    for(int i = 0; i < instances; i++) {
        colorsTransforms[i].v[0] = (float)colors[i].r / 255.0f;
        colorsTransforms[i].v[1] = (float)colors[i].g / 255.0f;
        colorsTransforms[i].v[2] = (float)colors[i].b / 255.0f;
        colorsTransforms[i].v[3] = (float)colors[i].a / 255.0f;
    }

    // Enable mesh VAO to attach new buffer
    rlEnableVertexArray(mesh.vaoId);

    // This could alternatively use a static VBO and either glMapBuffer() or glBufferSubData().
    // It isn't clear which would be reliably faster in all cases and on all platforms,
    // anecdotally glMapBuffer() seems very slow (syncs) while glBufferSubData() seems
    // no faster, since we're transferring all the transform matrices anyway
    instancesVboId = rlLoadVertexBuffer(instanceTransforms, (int)(instances * sizeof(float16)), false);

    // Instances transformation matrices are send to shader attribute location: SHADER_LOC_MATRIX_MODEL
    for(unsigned int i = 0; i < 4; i++) {
        rlEnableVertexAttribute(shader.locs[SHADER_LOC_MATRIX_MODEL] + i);
        rlSetVertexAttribute(shader.locs[SHADER_LOC_MATRIX_MODEL] + i, 4, RL_FLOAT, 0, sizeof(Matrix), (void *)(i * sizeof(Vector4)));
        rlSetVertexAttributeDivisor(shader.locs[SHADER_LOC_MATRIX_MODEL] + i, 1);
    }

    rlDisableVertexBuffer();
    rlDisableVertexArray();

    // Enable mesh VAO to attach new buffer
    rlEnableVertexArray(mesh.vaoId);
    colorsVboId = rlLoadVertexBuffer(colorsTransforms, (int)(instances * sizeof(float4)), true);

    // Colors are send to shader attribute location: SHADER_LOC_VERTEX_COLOR
    rlEnableVertexAttribute(shader.locs[SHADER_LOC_VERTEX_COLOR]);
    rlSetVertexAttribute(shader.locs[SHADER_LOC_VERTEX_COLOR], 4, RL_FLOAT, 0, sizeof(float4), 0);
    rlSetVertexAttributeDivisor(shader.locs[SHADER_LOC_VERTEX_COLOR], 1);

    rlDisableVertexBuffer();
    rlDisableVertexArray();

    // Accumulate internal matrix transform (push/pop) and view matrix
    // NOTE: In this case, model instance transformation must be computed in the shader
    matModelView = MatrixMultiply(rlGetMatrixTransform(), matView);

    // Upload model normal matrix (if locations available)
    if(shader.locs[SHADER_LOC_MATRIX_NORMAL] != -1)
        rlSetUniformMatrix(shader.locs[SHADER_LOC_MATRIX_NORMAL], MatrixTranspose(MatrixInvert(matModel)));

    int dgrid_loc = GetShaderLocation(shader, "dgrid");
    if(dgrid_loc != 1) {
        rlSetUniform(dgrid_loc, (void *)&grid_mask, RL_SHADER_UNIFORM_INT, 1);
    }

    rlEnableVertexArray(mesh.vaoId);
    if(mesh.indices != NULL)
        rlEnableVertexBufferElement(mesh.vboId[6]);

    // Calculate model-view-projection matrix (MVP)
    Matrix matModelViewProjection;
    matModelViewProjection = MatrixMultiply(matModelView, matProjection);

    // Send combined model-view-projection matrix to shader
    rlSetUniformMatrix(shader.locs[SHADER_LOC_MATRIX_MVP], matModelViewProjection);

    // Draw mesh instanced
    if(mesh.indices != NULL) {
        rlDrawVertexArrayElementsInstanced(0, mesh.triangleCount * 3, 0, instances);
    } else {
        rlDrawVertexArrayInstanced(0, mesh.vertexCount, instances);
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
