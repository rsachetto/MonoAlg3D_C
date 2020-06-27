set(SOURCE_FILES_CONFIG
        config_parser.h config_parser.cpp)

ADD_LIBRARY(config_parser STATIC ${SOURCE_FILES_CONFIG})

set_property(TARGET config_parser PROPERTY POSITION_INDEPENDENT_CODE ON)
