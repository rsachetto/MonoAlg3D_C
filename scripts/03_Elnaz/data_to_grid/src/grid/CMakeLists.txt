set(SOURCE_FILES_GRID
          grid.h grid.cpp
          hexaedron_cell.h hexaedron_cell.cpp)


ADD_LIBRARY(grid STATIC ${SOURCE_FILES_GRID})

set_property(TARGET grid PROPERTY POSITION_INDEPENDENT_CODE ON)
