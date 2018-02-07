cmake_minimum_required(VERSION 3.0.1)
project(MonoAlg3D)

IF(MSVC)
	set(GLUT_INCLUDE_DIR "C:/Users/sache/freeglut/include")
	SET(GLUT_glut_LIBRARY "C:/Users/sache/freeglut/lib/x64/freeglut.lib")
endif()

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
    ENDIF(CMAKE_C_COMPILER_ID)
ELSE()
    SET(CMAKE_C_STANDARD 99)
ENDIF()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )

IF(EXISTS "/etc/manjaro-release")
    set(CMAKE_C_COMPILER "/opt/cuda/bin/gcc")
ENDIF(EXISTS "/etc/manjaro-release")

IF(MSVC)
	SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
	ADD_SUBDIRECTORY(src/getopt)
	ADD_SUBDIRECTORY(src/dlfcn-win32)
	SET(OPT_DEPS getopt)
ELSE()
	SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing -Wall -Wno-unused-function")
ENDIF()

FIND_PACKAGE(OpenMP)
IF (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF(OPENMP_FOUND)

FIND_PACKAGE(CUDA)

IF (CUDA_FOUND AND MSVC AND "${CMAKE_SIZEOF_VOID_P}" EQUAL "4")
            message(WARNING "You have CUDA installed, but we can't use it unless you put visual studio in 64bit mode.")
            set(CUDA_FOUND 0)
ENDIF()

IF (CUDA_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I \"${CUDA_INCLUDE_DIRS}\" -DCOMPILE_CUDA")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I \"${CUDA_INCLUDE_DIRS}\" -DCOMPILE_CUDA")

    #windos hack
	IF (MSVC OR APPLE)
        set(CUDA_PROPAGATE_HOST_FLAGS OFF)
       # Grab all the -D flags from CMAKE_CXX_FLAGS so we can pass them
       # to nvcc.
       string(REGEX MATCHALL "-D[^ ]*" FLAGS_FOR_NVCC "${CMAKE_CXX_FLAGS}")
       IF(MSVC)
           list(APPEND CUDA_NVCC_FLAGS "${FLAGS_FOR_NVCC} -Xcompiler \"-MDd  ${OpenMP_CXX_FLAGS}\" ")
       ENDIF()
    ENDIF()

    ADD_SUBDIRECTORY(src/gpu_utils)
    SET(OPT_DEPS gpu_utils)
ENDIF(CUDA_FOUND)

FIND_PACKAGE(OpenGL)
FIND_PACKAGE(GLUT)
IF(OPENGL_FOUND AND GLUT_FOUND)
	IF(MSVC)
		SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I \"${GLUT_INCLUDE_DIR}\" -DCOMPILE_OPENGL")
	ELSE()
		SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I \"${OPENGL_INCLUDE_DIR}\" -I \"${GLUT_INCLUDE_DIR}\" -DCOMPILE_OPENGL")		
	ENDIF()

    ADD_SUBDIRECTORY(src/draw)
    SET(OPT_DEPS ${OPT_DEPS} draw)
ENDIF()

ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/alg)
ADD_SUBDIRECTORY(src/monodomain)
ADD_SUBDIRECTORY(src/hash)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/string)

ADD_EXECUTABLE(MonoAlg3D src/main.c)

ADD_DEPENDENCIES(MonoAlg3D utils alg solvers hashes ini_parser string ${OPT_DEPS})

IF(MSVC)
    TARGET_LINK_LIBRARIES(MonoAlg3D alg solvers utils hashes ini_parser string ${OPT_DEPS} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${CUDA_LIBRARIES} dl)
    SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )
ELSE()
    TARGET_LINK_LIBRARIES(MonoAlg3D alg solvers utils hashes ini_parser string ${OPT_DEPS} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${CUDA_LIBRARIES} dl m)
ENDIF()


ADD_SUBDIRECTORY(src/models_library)

SET(COMMON_SRC ../libraries_common/config_helpers.c ../libraries_common/config_helpers.h
        ../libraries_common/erros_helpers.c ../libraries_common/erros_helpers.h)

ADD_SUBDIRECTORY(src/stimuli_library)
ADD_SUBDIRECTORY(src/domains_library)
ADD_SUBDIRECTORY(src/extra_data_library)
###

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests_bin )

ADD_SUBDIRECTORY(src/tests)
