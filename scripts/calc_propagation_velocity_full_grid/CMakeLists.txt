cmake_minimum_required(VERSION 2.8)
 
PROJECT(CalcPropagationVelocity)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/tissue)

ADD_EXECUTABLE(CalcPropagationVelocity src/main.cpp )

ADD_DEPENDENCIES(CalcPropagationVelocity utils tissue)
 
TARGET_LINK_LIBRARIES(CalcPropagationVelocity ${VTK_LIBRARIES} utils tissue)
