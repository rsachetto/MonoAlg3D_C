LIB_STATIC_LIBS="config_helpers utils tinyexpr"

COMPILE_SHARED_LIB "default_extra_data" "extra_data.c" "" "${LIB_STATIC_LIBS}"
COMPILE_SHARED_LIB "mixed_models_extra_data" "extra_mixed_models.c" "" "${LIB_STATIC_LIBS}"