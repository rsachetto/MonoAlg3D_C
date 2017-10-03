cmake_minimum_required(VERSION 3.8)
project(MonoAlg3D)

set(CMAKE_C_STANDARD 99)

if(EXISTS "/etc/manjaro-release")
    set(CMAKE_C_COMPILER "/opt/cuda/bin/gcc")
endif(EXISTS "/etc/manjaro-release")

find_package (Threads) #we look for a thread package

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(CUDA)

SET(CUDA_NVCC_FLAGS "-Xcompiler -fopenmp" CACHE STRING "nvcc flags" FORCE)
#HW capability 1.3
#SET (CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" FORCE)

set(GPU_FILE src/gpu/gpu_common.cu src/gpu/gpu_common.h)

SET(LIB_TYPE STATIC) #set the lib type
CUDA_ADD_LIBRARY(CudaLib ${LIB_TYPE} ${GPU_FILE})

set(SOURCE_FILES src/main.c
        src/grid/cell/cell.h src/grid/cell/cell.c
        src/grid/grid.c src/grid/grid.h
        src/solvers/monodomain_solver.h
        src/grid/grid_refinement.c src/grid/grid_derefinement.c
        src/grid/cell/cell_derefinement.c src/grid/cell/cell_refinement.c
        src/utils/opts.h src/grid/grid_domain.c src/utils/sort.c src/utils/utils.h src/utils/search.c
        src/utils/point_hash.c src/utils/point_hash.h
        src/solvers/ode_solver.c src/solvers/ode_solver.h
        src/solvers/monodomain_solver.h src/solvers/monodomain_solver.c
        src/utils/output_utils.h src/utils/output_utils.c
        src/solvers/constants.h
        src/utils/stop_watch.h src/utils/stop_watch.c)

add_executable(MonoAlg3D ${SOURCE_FILES})

ADD_DEPENDENCIES(MonoAlg3D CudaLib)
target_link_libraries(MonoAlg3D CudaLib ${CMAKE_THREAD_LIBS_INIT} m)