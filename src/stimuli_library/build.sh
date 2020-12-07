CHECK_CUSTOM_FILE "custom_stimuli_functions.c"

COMPILE_SHARED_LIB "default_stimuli" "stimuli.c ${CUSTOM_FILE}" "" "alg config_helpers utils tinyexpr" "m"
