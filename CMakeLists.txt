cmake_minimum_required(VERSION 3.0.1)
project(MonoAlg3D)

macro(USE_C99)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
            SET(CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
        endif(CMAKE_C_COMPILER_ID)
    else()
        SET(CMAKE_C_STANDARD 99)
    endif()
endmacro(USE_C99)

USE_C99()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )

if(EXISTS "/etc/manjaro-release")
    set(CMAKE_C_COMPILER "/opt/cuda/bin/gcc")
endif(EXISTS "/etc/manjaro-release")

SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing -Wall -Wno-unused-function")

FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

FIND_PACKAGE(CUDA)
if (CUDA_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I ${CUDA_INCLUDE_DIRS} -DCOMPILE_CUDA")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CUDA_INCLUDE_DIRS} -DCOMPILE_CUDA")
endif()

FIND_PACKAGE(OpenGL)
FIND_PACKAGE(GLUT)
if(OPENGL_FOUND AND GLUT_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I ${OPENGL_INCLUDE_DIR} -I ${GLUT_INCLUDE_DIR} -DCOMPILE_OPENGL")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${OPENGL_INCLUDE_DIR} -I ${GLUT_INCLUDE_DIR} -DCOMPILE_OPENGL")
    ADD_SUBDIRECTORY(src/draw)
endif()

ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/gpu_utils)
ADD_SUBDIRECTORY(src/alg)
ADD_SUBDIRECTORY(src/monodomain)
ADD_SUBDIRECTORY(src/hash)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/string)

SET(SOURCE_FILES src/main.c)

ADD_EXECUTABLE(MonoAlg3D ${SOURCE_FILES})

ADD_DEPENDENCIES(MonoAlg3D utils alg solvers gpu_utils hashes ini_parser string)

if(OPENGL_FOUND AND GLUT_FOUND)
    TARGET_LINK_LIBRARIES(MonoAlg3D alg solvers utils gpu_utils hashes ini_parser string draw ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${CUDA_LIBRARIES} dl m)
else()
    TARGET_LINK_LIBRARIES(MonoAlg3D alg solvers utils gpu_utils hashes ini_parser string  ${CUDA_LIBRARIES} dl m)
endif()


ADD_SUBDIRECTORY(src/models_library)
ADD_SUBDIRECTORY(src/stimuli_library)
ADD_SUBDIRECTORY(src/domains_library)
ADD_SUBDIRECTORY(src/extra_data_library)

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests_bin )

ADD_SUBDIRECTORY(src/tests)
