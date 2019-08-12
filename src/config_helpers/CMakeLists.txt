set(SOURCE_FILES_CONFIG config_helpers.h config_helpers.c)

ADD_LIBRARY(config_helpers STATIC ${SOURCE_FILES_CONFIG})
set_property(TARGET config_helpers PROPERTY POSITION_INDEPENDENT_CODE ON)