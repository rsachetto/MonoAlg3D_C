cmake_minimum_required(VERSION 3.0.1)
project(MonoAlg3D)

macro(use_c99)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
        set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    endif ()
    else ()
        set (CMAKE_C_STANDARD 99)
    endif ()
endmacro(use_c99)


#set(CMAKE_C_STANDARD 99)
use_c99()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )

if(EXISTS "/etc/manjaro-release")
    set(CMAKE_C_COMPILER "/opt/cuda/bin/gcc")
endif(EXISTS "/etc/manjaro-release")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing -Wall -Wno-unused-function")

find_package (Threads) #we look for a thread package

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

FIND_PACKAGE(CUDA)
if (CUDA_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I ${CUDA_INCLUDE_DIRS} -DCOMPILE_CUDA")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CUDA_INCLUDE_DIRS} -DCOMPILE_CUDA")
endif()

ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/gpu_utils)
ADD_SUBDIRECTORY(src/alg)
ADD_SUBDIRECTORY(src/main)
ADD_SUBDIRECTORY(src/hash)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/string)

set(SOURCE_FILES src/main.c)

add_executable(MonoAlg3D ${SOURCE_FILES})

ADD_DEPENDENCIES(MonoAlg3D utils alg solvers gpu_utils hashes ini_parser string)

target_link_libraries(MonoAlg3D alg solvers utils gpu_utils hashes ini_parser string ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES} dl m)

ADD_SUBDIRECTORY(src/models_library)
ADD_SUBDIRECTORY(src/stimuli_library)
ADD_SUBDIRECTORY(src/domains_library)
ADD_SUBDIRECTORY(src/extra_data_library)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests_bin )

ADD_SUBDIRECTORY(src/tests)
