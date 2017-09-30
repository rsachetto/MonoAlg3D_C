cmake_minimum_required(VERSION 3.8)
project(MonoAlg3D)

set(CMAKE_C_STANDARD 99)

find_package (Threads)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES src/main.c
        src/cell.h src/cell.c
        src/grid.c src/grid.h
        src/solver.h
        src/grid_refinement.c src/grid_derefinement.c
        src/cell_derefinement.c src/cell_refinement.c
        src/opts.h)

add_executable(MonoAlg3D ${SOURCE_FILES})
target_link_libraries(MonoAlg3D ${CMAKE_THREAD_LIBS_INIT})