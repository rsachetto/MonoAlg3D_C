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
        src/grid/cell/cell.h src/grid/cell/cell.c
        src/grid/grid.c src/grid/grid.h
        src/solver.h
        src/grid/grid_refinement.c src/grid/grid_derefinement.c
        src/grid/cell/cell_derefinement.c src/grid/cell/cell_refinement.c
        src/opts.h src/grid/grid_domain.c src/utils/sort.c src/utils/utils.h src/utils/search.c src/utils/point_hash.c src/utils/point_hash.h)

add_executable(MonoAlg3D ${SOURCE_FILES})
target_link_libraries(MonoAlg3D ${CMAKE_THREAD_LIBS_INIT} m)