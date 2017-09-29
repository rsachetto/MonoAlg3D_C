cmake_minimum_required(VERSION 3.8)
project(MonoAlg3D)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES src/main.c src/cell.h src/cell.c src/grid.c src/grid.h src/solver.h)
add_executable(MonoAlg3D ${SOURCE_FILES})