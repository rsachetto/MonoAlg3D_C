LIB_STATIC_LIBS="config_helpers utils tinyexpr"

CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_extra_data" "extra_data.c helper_functions.c ${CUSTOM_FILE}" "helper_functions.h" "${LIB_STATIC_LIBS}"
