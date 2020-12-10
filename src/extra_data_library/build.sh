LIB_STATIC_LIBS="config_helpers utils tinyexpr"

CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_extra_data" "extra_data.c ${CUSTOM_FILE}" "" "${LIB_STATIC_LIBS}"
