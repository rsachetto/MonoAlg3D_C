LIB_STATIC_DEPS="utils alg config_helpers tinyexpr"

COMPILE_SHARED_LIB "default_matrix_assembly" "matrix_assembly.c" "" "${LIB_STATIC_DEPS}"

if [ -n "$COMPILE_WITH_DDM" ]; then
    COMPILE_SHARED_LIB "ddm_matrix_assembly" "ddm_matrix_assembly.c" "" "${LIB_STATIC_DEPS}"
fi

COMPILE_SHARED_LIB "purkinje_matrix_assembly" "purkinje_matrix_assembly.c" "" "${LIB_STATIC_DEPS}"

COMPILE_SHARED_LIB "purkinje_coupled_matrix_assembly" "purkinje_coupled_matrix_assembly.c" "" "${LIB_STATIC_DEPS}"
