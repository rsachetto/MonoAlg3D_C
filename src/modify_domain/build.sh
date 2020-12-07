CHECK_CUSTOM_FILE "custom_modify_domain_functions.c"

COMPILE_SHARED_LIB "default_modify_domain" "modify_domain.c ${CUSTOM_FILE}" "" "tinyexpr"
