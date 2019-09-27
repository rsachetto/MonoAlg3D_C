set(SOURCE_FILES_TISSUE
        tissue.h tissue.cpp
        cell.h cell.cpp)

ADD_LIBRARY(tissue STATIC ${SOURCE_FILES_TISSUE})

set_property(TARGET tissue PROPERTY POSITION_INDEPENDENT_CODE ON)
