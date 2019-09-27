cmake_minimum_required(VERSION 2.8)
 
PROJECT(PurkinjeReader)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

ADD_SUBDIRECTORY(graph)
ADD_SUBDIRECTORY(purkinje)

add_executable(PurkinjeReader main.cpp )
ADD_DEPENDENCIES(PurkinjeReader graph purkinje)
 
target_link_libraries(PurkinjeReader graph purkinje ${VTK_LIBRARIES})
