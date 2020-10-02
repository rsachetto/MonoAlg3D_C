SAVE_STATIC_DEPS="vtk_utils config_helpers utils sds alg tinyexpr miniz gpu_utils"

if [ -n "$CUDA_FOUND" ]; then
  SAVE_STATIC_DEPS="$SAVE_STATIC_DEPS gpu_utils"
  SAVE_DYNAMIC_DEPS="cudart"
fi


COMPILE_SHARED_LIB "default_save_mesh" "save_mesh.c" "" "$SAVE_STATIC_DEPS" "$SAVE_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH"

COMPILE_SHARED_LIB "default_save_mesh_purkinje" "save_mesh_purkinje.c save_mesh_helper.c" "save_mesh_helper.h" "$SAVE_STATIC_DEPS" "$SAVE_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH"
