CHECK_CUSTOM_FILE "custom_purkinje_functions.c"

COMPILE_SHARED_LIB "default_purkinje" "purkinje.c ${CUSTOM_FILE} purkinje_helpers.c" "purkinje_helpers.h" "alg config_helpers utils solvers graph tinyexpr"
