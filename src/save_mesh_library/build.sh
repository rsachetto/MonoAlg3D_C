SAVE_STATIC_DEPS="vtk_utils ensight_utils graph config_helpers utils sds alg tinyexpr miniz"

if [ -n "$CUDA_FOUND" ]; then
  SAVE_STATIC_DEPS="$SAVE_STATIC_DEPS"
  SAVE_DYNAMIC_DEPS="cudart"
fi

CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_save_mesh" "save_mesh_helper.c save_mesh.c ${CUSTOM_FILE}" "save_mesh_helper.h" "$SAVE_STATIC_DEPS" "$SAVE_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH"

COMPILE_SHARED_LIB "default_save_mesh_purkinje" "save_mesh_helper.c save_mesh_purkinje.c" "save_mesh_helper.h" "$SAVE_STATIC_DEPS" "$SAVE_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH"
