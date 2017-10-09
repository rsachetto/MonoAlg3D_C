cmake_minimum_required(VERSION 3.8)
project(MonoAlg3D)

set(CMAKE_C_STANDARD 99)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/model_lib )

if(EXISTS "/etc/manjaro-release")
    set(CMAKE_C_COMPILER "/opt/cuda/bin/gcc")
endif(EXISTS "/etc/manjaro-release")

find_package (Threads) #we look for a thread package

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

FIND_PACKAGE(CUDA)

ADD_SUBDIRECTORY(src/models)
ADD_SUBDIRECTORY(src/gpu)
ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/grid)
#Maibe solver should not be a independent library as it depends on utils and gpu
ADD_SUBDIRECTORY(src/solvers)

set(SOURCE_FILES src/main.c)

add_executable(MonoAlg3D ${SOURCE_FILES})

ADD_DEPENDENCIES(MonoAlg3D CudaLib utils alg solvers)

target_link_libraries(MonoAlg3D alg solvers utils CudaLib ${CMAKE_THREAD_LIBS_INIT} m)

ADD_SUBDIRECTORY(src/tests)
