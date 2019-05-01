cmake_minimum_required(VERSION 2.8)
 
PROJECT(CalcPropagationVelocity)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

add_executable(CalcPropagationVelocity main.cpp )
 
target_link_libraries(CalcPropagationVelocity ${VTK_LIBRARIES})
