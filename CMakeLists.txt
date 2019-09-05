IF(EXISTS "/etc/manjaro-release")
    set(CMAKE_C_COMPILER "/opt/cuda/bin/gcc")
    set(CMAKE_CXX_COMPILER "/opt/cuda/bin/g++")
ENDIF(EXISTS "/etc/manjaro-release")

#IF(EXISTS "/etc/debconf.conf")
#    set(CMAKE_C_COMPILER "/usr/bin/gcc-6")
#    set(CMAKE_CXX_COMPILER "/usr/bin/g++-6")
#ENDIF(EXISTS "/etc/debconf.conf")

cmake_minimum_required(VERSION 3.0.1)
project(MonoAlg3D)

set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

if (CMAKE_VERSION VERSION_LESS "3.9.0")
    # Set only because there are old target_link_libraries() statements in the
    # FindCUDA.cmake file that comes with CMake that error out if the new behavior
    # is used.  In newer versions of CMake we can instead set CUDA_LINK_LIBRARIES_KEYWORD which fixes this issue.
    cmake_policy(SET CMP0023 OLD)
else()
    set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
endif()

IF (CMAKE_VERSION VERSION_LESS "3.1")
    IF (CMAKE_C_COMPILER_ID STREQUAL "GNU")
        SET(CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    ENDIF(CMAKE_C_COMPILER_ID STREQUAL "GNU")
ELSE()
    SET(CMAKE_C_STANDARD 99)
ENDIF()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )

SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing -Wall -Wno-unused-function -Wno-char-subscripts")

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -DDEBUG_INFO" )
ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")

FIND_PACKAGE(OpenMP)
IF (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF(OPENMP_FOUND)

FIND_PACKAGE(CUDA)

IF (CUDA_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I \"${CUDA_INCLUDE_DIRS}\" -DCOMPILE_CUDA")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I \"${CUDA_INCLUDE_DIRS}\" -DCOMPILE_CUDA")

	IF (APPLE)
        set(CUDA_PROPAGATE_HOST_FLAGS OFF)
       # Grab all the -D flags from CMAKE_CXX_FLAGS so we can pass them
       # to nvcc.
       string(REGEX MATCHALL "-D[^ ]*" FLAGS_FOR_NVCC "${CMAKE_CXX_FLAGS}")
    ENDIF()
ENDIF(CUDA_FOUND)

SET(OpenGL_GL_PREFERENCE "GLVND")
FIND_PACKAGE(OpenGL)

if (OPENGL_FOUND AND OPENMP_FOUND)
    set (COMPILE_GUI TRUE)
else ()
    set (COMPILE_GUI FALSE)
endif ()

IF(COMPILE_GUI)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DCOMPILE_OPENGL")
ENDIF(COMPILE_GUI)

ADD_SUBDIRECTORY(src/config_helpers)
ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/alg)
ADD_SUBDIRECTORY(src/monodomain)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/string)
ADD_SUBDIRECTORY(src/config)
ADD_SUBDIRECTORY(src/graph)
ADD_SUBDIRECTORY(src/vtk_utils)

IF(COMPILE_GUI)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(BUILD_GAMES    OFF CACHE BOOL "" FORCE)
    ADD_SUBDIRECTORY(src/raylib/src)
    ADD_SUBDIRECTORY(src/draw)
    SET(OPT_DEPS ${OPT_DEPS} draw raylib)
ENDIF(COMPILE_GUI)

IF (CUDA_FOUND)
    ADD_SUBDIRECTORY(src/gpu_utils)
    SET(OPT_DEPS ${OPT_DEPS} gpu_utils)
ENDIF()

ADD_EXECUTABLE(MonoAlg3D src/main_simulator.c)
ADD_DEPENDENCIES(MonoAlg3D utils alg solvers ini_parser string graph config ${OPT_DEPS})

FIND_PACKAGE(MPI)
IF(MPI_C_FOUND)
    INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
    ADD_EXECUTABLE(MonoAlg3D_batch src/main_batch.c)
    ADD_DEPENDENCIES(MonoAlg3D_batch utils alg solvers ini_parser string graph config ${OPT_DEPS})
    TARGET_LINK_LIBRARIES(MonoAlg3D_batch alg solvers utils ini_parser string graph config ${OPT_DEPS} ${CUDA_LIBRARIES} ${MPI_C_LIBRARIES}  dl m)
ENDIF(MPI_C_FOUND)

TARGET_LINK_LIBRARIES(MonoAlg3D alg solvers utils ini_parser string graph config ${OPT_DEPS} ${CUDA_LIBRARIES} dl m)

ADD_SUBDIRECTORY(src/models_library)
ADD_SUBDIRECTORY(src/stimuli_library)
ADD_SUBDIRECTORY(src/domains_library)
ADD_SUBDIRECTORY(src/purkinje_library)
ADD_SUBDIRECTORY(src/extra_data_library)
ADD_SUBDIRECTORY(src/matrix_assembly_library)
ADD_SUBDIRECTORY(src/linear_system_solver_library)
ADD_SUBDIRECTORY(src/save_mesh_library)
ADD_SUBDIRECTORY(src/save_state_library)
ADD_SUBDIRECTORY(src/restore_state_library)
ADD_SUBDIRECTORY(src/update_monodomain_library)
ADD_SUBDIRECTORY(src/xml_parser)

#BUILDING VISUALIZER
IF(COMPILE_GUI)
    ADD_EXECUTABLE(MonoAlg3D_visualizer src/main_visualizer.c)
    ADD_DEPENDENCIES(MonoAlg3D_visualizer utils config vtk_utils draw)
    TARGET_LINK_LIBRARIES(MonoAlg3D_visualizer utils config vtk_utils draw raylib)

#    ADD_EXECUTABLE(MonoAlg3D_activation_map_analyzer src/main_activation_map_analyzer.c)
#    ADD_DEPENDENCIES(MonoAlg3D_activation_map_analyzer alg utils config vtk_utils)
#    TARGET_LINK_LIBRARIES(MonoAlg3D_activation_map_analyzer alg utils config vtk_utils)
ENDIF(COMPILE_GUI)

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests_bin )

ADD_SUBDIRECTORY(src/tests)