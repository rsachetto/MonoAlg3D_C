cmake_minimum_required(VERSION 2.8)
 
PROJECT(ScenariosConfiguration)
 

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

add_executable(ScenariosConfiguration main.cpp )
 
target_link_libraries(ScenariosConfiguration)
