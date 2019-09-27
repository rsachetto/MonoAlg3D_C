set(SOURCE_FILES_UTILS
        utils.h utils.cpp
        file_utils.h file_utils.cpp)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})

set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)
