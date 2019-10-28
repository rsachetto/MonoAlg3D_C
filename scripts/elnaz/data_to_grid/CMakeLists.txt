cmake_minimum_required(VERSION 2.8)

PROJECT(DataToGrid)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/config)
ADD_SUBDIRECTORY(src/grid)
#ADD_SUBDIRECTORY(src/tissue)

ADD_EXECUTABLE(DataToGrid src/main.cpp )

ADD_DEPENDENCIES(DataToGrid utils config_parser grid)

TARGET_LINK_LIBRARIES(DataToGrid ${VTK_LIBRARIES} utils config_parser grid)
