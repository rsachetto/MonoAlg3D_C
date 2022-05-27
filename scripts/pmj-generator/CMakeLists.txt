cmake_minimum_required(VERSION 3.0.1)
 
PROJECT(PMJGenerator)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

ADD_SUBDIRECTORY(src/graph)
ADD_SUBDIRECTORY(src/purkinje)

add_executable(PMJGenerator src/main.cpp )
ADD_DEPENDENCIES(PMJGenerator graph purkinje)
 
target_link_libraries(PMJGenerator ${VTK_LIBRARIES} graph purkinje)
