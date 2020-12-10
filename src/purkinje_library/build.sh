CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_purkinje" "purkinje.c ${CUSTOM_FILE} purkinje_helpers.c" "purkinje_helpers.h" "alg config_helpers utils solvers graph tinyexpr"
