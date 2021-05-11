cmake_minimum_required(VERSION 2.8)
 
PROJECT(tuneCV)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

add_executable(tuneCV main.cpp )

target_link_libraries(tuneCV ${VTK_LIBRARIES})
