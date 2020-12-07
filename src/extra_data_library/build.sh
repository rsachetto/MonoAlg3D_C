LIB_STATIC_LIBS="config_helpers utils tinyexpr"

CHECK_CUSTOM_FILE "custom_extra_data_functions.c"

COMPILE_SHARED_LIB "default_extra_data" "extra_data.c ${CUSTOM_FILE}" "" "${LIB_STATIC_LIBS}"
