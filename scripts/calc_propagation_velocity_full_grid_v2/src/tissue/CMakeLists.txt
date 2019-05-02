set(SOURCE_FILES_TISSUE
        tissue.h tissue.cpp
        cell.h cell.cpp
        control_volume.h control_volume.cpp
        cell_data.h cell_data.cpp)

ADD_LIBRARY(tissue STATIC ${SOURCE_FILES_TISSUE})

set_property(TARGET tissue PROPERTY POSITION_INDEPENDENT_CODE ON)
