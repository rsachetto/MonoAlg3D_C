CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_calc_ecg" "ecg.c ${CUSTOM_FILE}" "" "config_helpers vtk_utils utils alg sds tinyexpr" "m"
