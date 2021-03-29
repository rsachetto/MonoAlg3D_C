CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_domains" "domain.c ${CUSTOM_FILE} domain_helpers.c" "domain_helpers.h" "config_helpers vtk_utils utils alg sds tinyexpr" "m"
