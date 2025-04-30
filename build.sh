#!/bin/bash

PRINT_USAGE () {
    echo "Usage $0 [flags] [modules]" >&2;
    echo "Valid modules: all, gui, simulator or batch (default is all)" >&2;
    echo "Valid flags:" >&2;
    echo "-f  - force recompilation" >&2;
    echo "-l  - write build log on compile_commands.json" >&2;
    echo "-q  - quiet compilation. Only errors and warnings will be outputed" >&2;
    echo "-r  - build release version (Default)" >&2;
    echo "-d  - build debug version" >&2;
    exit 1
}

OPTIONS_FILE=./bsbash/build_functions.sh
FUNCTIONS_FILE=./bsbash/find_functions.sh

if [ -f "$OPTIONS_FILE" ]; then
    # shellcheck disable=SC1090
    source $OPTIONS_FILE
else
    echo "$OPTIONS_FILE not found, aborting compilation"
    exit 1
fi

if [ -f "$FUNCTIONS_FILE" ]; then
    # shellcheck disable=SC1090
    source $FUNCTIONS_FILE
fi

CHECK_COLOR_SUPPORT

###########User code#####################
COMPILE_GUI=''
COMPILE_MPI=''
COMPILE_CONVERTER=''
COMPILE_FIBER_CONVERTER=''
COMPILE_SIMULATOR=''
COMPILE_POSTPROCESSOR=''
COMPILE_WITH_DDM=''
DISABLE_CUDA=''
USE_SYCL=''

GET_BUILD_OPTIONS "$@"

#Comment this to compile with single precision floats on the GPU
C_FLAGS="-DCELL_MODEL_REAL_DOUBLE"

if [ "$BUILD_TYPE" == "release" ]; then
    C_FLAGS="$C_FLAGS -O3"
elif [ "$BUILD_TYPE" == "debug" ]; then
    #C_FLAGS="$C_FLAGS -g -DDEBUG_INFO -fprofile-arcs -ftest-coverage"
    #C_FLAGS="$C_FLAGS -g -DDEBUG_INFO -fsanitize=address"
    C_FLAGS="$C_FLAGS -g -DDEBUG_INFO"
else
    PRINT_ERROR "$BUILD_TYPE is not a valid BUILD_TYPE."
    PRINT_ERROR "Valid BUILD_TYPE options are: release, debug (-r or -d options)"
    exit 1
fi

for i in "${BUILD_ARGS[@]}"; do
    case $i in
        clean)
            echo "Cleaning $BUILD_TYPE"
            CLEAN_PROJECT "$BUILD_TYPE"
            rm -fr shared_libs/
            cd src/3dparty/raylib/src || exit 1;
            make clean
            cd - || exit 1;
            exit 0
            ;;
        ddm)
            C_FLAGS="$C_FLAGS -DENABLE_DDM"
            COMPILE_WITH_DDM='y'
            COMPILE_GUI='y'
            COMPILE_MPI='y'
            COMPILE_SIMULATOR='y'
            COMPILE_CONVERTER='y'
            COMPILE_FIBER_CONVERTER='y'
            COMPILE_POSTPROCESSOR='y'
            ;;
        all)
            COMPILE_GUI='y'
            COMPILE_MPI='y'
            COMPILE_SIMULATOR='y'
            COMPILE_CONVERTER='y'
            COMPILE_FIBER_CONVERTER='y'
            COMPILE_POSTPROCESSOR='y'
            COMPILE_EXPAND='y'
            COMPILE_CLIP='y'
            COMPILE_EIKONAL='y'
            ;;
        simulator)
            COMPILE_SIMULATOR='y'
            ;;
        gui)
            COMPILE_GUI='y'
            ;;
        batch)
            COMPILE_MPI='y'
            ;;
        converter)
            COMPILE_CONVERTER='y'
            ;;
        only_cpu)
            COMPILE_SIMULATOR='y'
            DISABLE_CUDA='y'
            ;;
        sycl)
            COMPILE_SIMULATOR='y'
	    COMPILE_MPI='y'
            DISABLE_CUDA='y'
            USE_SYCL='y'
            ;;
        *)
            echo "Invalid option $i. Aborting!"
            exit 1
            ;;
    esac
done

DEFAULT_C_FLAGS="-fopenmp -std=gnu99 -fno-strict-aliasing  -Wall -Wno-stringop-truncation -Wno-unused-function -Wno-char-subscripts -Wno-unused-result -Wno-switch -Werror=implicit-function-declaration"

# Update the SYCL compiler and Intel flags
if [ -n "$USE_SYCL" ]; then
    C_COMPILER="icx"
    CXX_COMPILER="icpx"
    DEFAULT_C_FLAGS="-qopenmp -fno-strict-aliasing  -Wall -Wno-unused-function -Wno-unused-result -Wno-switch -Werror=implicit-function-declaration"
fi

RUNTIME_OUTPUT_DIRECTORY="$ROOT_DIR/bin"
LIBRARY_OUTPUT_DIRECTORY="$ROOT_DIR/shared_libs"

C_FLAGS="$C_FLAGS $DEFAULT_C_FLAGS"

GET_LINUX_VERSION

echo -e "${INFO}Linux version: ${OS}-${VER}"

if [ -n "$COMPILE_SIMULATOR" ] || [ -n "$COMPILE_MPI" ]; then

    if [ -z "$DISABLE_CUDA" ]; then
       FIND_CUDA
    fi

    if [ -n "$USE_SYCL" ]; then
       FIND_SYCL
    fi

    echo -e "${INFO}C compiler:${NC} $C_COMPILER"
    echo -e "${INFO}C++ compiler:${NC} $CXX_COMPILER"

    if [ -n "$CUDA_FOUND" ]; then
        echo -e "${INFO}CUDA compiler:${NC} $NVCC"
        echo -e "${INFO}CUDA libraries path:${NC} $CUDA_LIBRARY_PATH"
        echo -e "${INFO}CUDA include path:${NC} $CUDA_INCLUDE_PATH"


        C_FLAGS="${C_FLAGS} -DCOMPILE_CUDA -I${CUDA_INCLUDE_PATH}"

    fi

    if [ -n "$SYCL_FOUND" ]; then
        echo -e "${INFO}SYCL compiler:${NC} ${ICPX}"

        C_FLAGS="${C_FLAGS} -DCOMPILE_SYCL -I${SYCL_INCLUDE_PATH}"

    fi
fi

if [ -n "$COMPILE_GUI" ]; then
        C_FLAGS="${C_FLAGS} -DCOMPILE_GUI"
fi

if [ -n "$CUDA_FOUND" ]; then
    if [ "$OS" == "Manjaro Linux" ]; then
        C_COMPILER=/opt/cuda/bin/gcc
        CXX_COMPILER=/opt/cuda/bin/g++
    elif [ "$OS" == "Ubuntu" ]; then
        if [ "$VER" == "20.10" ]; then
            C_COMPILER=gcc-9
            CXX_COMPILER=g++-9
        else
            C_COMPILER=gcc
            CXX_COMPILER=g++
        fi
    elif [ "$OS" == "Fedora" ]; then
        C_COMPILER=/usr/local/cuda/bin/gcc
        CXX_COMPILER=/usr/local/cuda/bin/g++
    fi

    FIND_AMGX
    if [ -n "$AMGX_FOUND" ]; then
        C_FLAGS="${C_FLAGS} -DAMGX -I${AMGX_INCLUDE_PATH}"
    fi

fi

echo -e "${INFO}C FLAGS:${NC} $C_FLAGS"

ADD_SUBDIRECTORY "src/3dparty/sds"
ADD_SUBDIRECTORY "src/config_helpers"
ADD_SUBDIRECTORY "src/utils"
ADD_SUBDIRECTORY "src/alg"
ADD_SUBDIRECTORY "src/monodomain"
ADD_SUBDIRECTORY "src/ode_solver"
ADD_SUBDIRECTORY "src/3dparty/ini_parser"
ADD_SUBDIRECTORY "src/config"
ADD_SUBDIRECTORY "src/graph"
ADD_SUBDIRECTORY "src/3dparty/xml_parser"
ADD_SUBDIRECTORY "src/3dparty/tinyexpr"
ADD_SUBDIRECTORY "src/3dparty/miniz"
ADD_SUBDIRECTORY "src/vtk_utils"
ADD_SUBDIRECTORY "src/ensight_utils"
ADD_SUBDIRECTORY "src/eikonal/"


#DINAMIC DEPS
ADD_SUBDIRECTORY "src/logger"

if [ -n "$COMPILE_GUI" ]; then
    ADD_SUBDIRECTORY "src/3dparty/raylib/src"
    ADD_SUBDIRECTORY "src/3dparty/tinyfiledialogs"
    ADD_SUBDIRECTORY "src/gui"
    OPT_DEPS_GUI="gui raylib tinyfiledialogs"
fi

if [ -n "$CUDA_FOUND" ] || [ -n "$SYCL_FOUND" ]; then
    ADD_SUBDIRECTORY "src/gpu_utils"
    OPT_DEPS_GPU=gpu_utils
fi

SRC_FILES="src/main_simulator.c"
HDR_FILES=""

STATIC_DEPS="monodomain ode_solver ini_parser config tinyexpr ${OPT_DEPS_GUI} config_helpers ensight_utils vtk_utils yxml alg graph utils sds miniz"

DYNAMIC_DEPS="dl m stdc++ ${OPT_DEPS_GPU} $CUDA_LIBRARIES"

if [ -n "$AMGX_FOUND" ]; then
    DYNAMIC_DEPS="$DYNAMIC_DEPS $AMGX_LIBRARIES"
fi


if [ -n "$COMPILE_GUI" ]; then
    DYNAMIC_DEPS="$DYNAMIC_DEPS OpenGL GLX GLU pthread X11 rt"
fi

DYNAMIC_DEPS="$DYNAMIC_DEPS logger"

if [ -n "$COMPILE_SIMULATOR" ] || [ -n "$COMPILE_BATCH" ]; then
    ADD_SUBDIRECTORY "src/models_library"
    ADD_SUBDIRECTORY "src/stimuli_library"
    ADD_SUBDIRECTORY "src/domains_library"
    ADD_SUBDIRECTORY "src/purkinje_library"
    ADD_SUBDIRECTORY "src/extra_data_library"
    ADD_SUBDIRECTORY "src/matrix_assembly_library"
    ADD_SUBDIRECTORY "src/linear_system_solver_library"
    ADD_SUBDIRECTORY "src/save_mesh_library"
    ADD_SUBDIRECTORY "src/save_state_library"
    ADD_SUBDIRECTORY "src/restore_state_library"
    ADD_SUBDIRECTORY "src/update_monodomain_library"
    ADD_SUBDIRECTORY "src/modify_domain"
    ADD_SUBDIRECTORY "src/ecg_library/"
fi

#COMPILE THE EXECUTABLES NOW
EXECUTABLES_LIBRARY_PATH="$CUDA_LIBRARY_PATH $CUDA_MATH_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

if [ -n "$AMGX_FOUND" ]; then
	EXECUTABLES_LIBRARY_PATH="$EXECUTABLES_LIBRARY_PATH $AMGX_LIBRARY_PATH"
fi



if [ -n "$COMPILE_SIMULATOR" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D" "$SRC_FILES" "$HDR_FILES" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH" "" "$USE_SYCL"
fi

if [ -n "$COMPILE_MPI" ]; then

  FIND_MPI

  if [ -n "$MPI_FOUND" ]; then
      SRC_FILES="src/main_batch.c"
      HDR_FILES=""
      DYNAMIC_DEPS="$DYNAMIC_DEPS $MPI_LIBRARIES"
      EXTRA_LIB_PATH="$EXTRA_LIB_PATH $MPI_LIBRARY_PATH"
    
      if [ -z "$MPI_INCLUDE_PATH" ]; then
	 COMPILE_EXECUTABLE "MonoAlg3D_batch" "$SRC_FILES" "$HDR_FILES" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH" "$INCLUDE_P" "$USE_SYCL"
      else
	 COMPILE_EXECUTABLE "MonoAlg3D_batch" "$SRC_FILES" "$HDR_FILES" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH" "$INCLUDE_P -I$MPI_INCLUDE_PATH" "$USE_SYCL"
      fi

  fi

fi

if [ -n "$COMPILE_GUI" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_visualizer" "src/main_visualizer.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH"
fi

if [ -n "$COMPILE_CONVERTER" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_converter" "src/main_converter.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH"
fi

if [ -n "$COMPILE_FIBER_CONVERTER" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_fiber_converter" "src/main_fiber_converter.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH"
fi

if [ -n "$COMPILE_POSTPROCESSOR" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_postprocessor" "src/main_postprocessor.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH"
    ADD_SUBDIRECTORY "src/postprocessing_library/"
fi

if [ -n "$COMPILE_EXPAND" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_expand_scar" "src/main_expand_scar.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH"
fi

if [ -n "$COMPILE_CLIP" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_clip_mesh" "src/main_clip_mesh.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH"
fi

if [ -n "$COMPILE_EIKONAL" ]; then
    COMPILE_EXECUTABLE "MonoAlg3D_eikonal_solver" "src/main_eikonal.c" "" "$STATIC_DEPS" "$DYNAMIC_DEPS eikonal_solver" "$EXECUTABLES_LIBRARY_PATH $EXTRA_LIB_PATH"
fi


FIND_CRITERION

if [ -n "$CRITERION_FOUND" ]; then
    # shellcheck disable=SC2034
    RUNTIME_OUTPUT_DIRECTORY=$ROOT_DIR/tests_bin
    ADD_SUBDIRECTORY "src/tests"
fi

if [ -n "$WAIT_ENTER" ]; then
    read -r -p "Press any key to continue... " -n1 -s
fi
