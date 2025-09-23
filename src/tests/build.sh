
if [ -n "$CUDA_FOUND" ] || [ -n "$SYCL_FOUND" ]; then
    TEST_OPT_DEPS=gpu_utils
fi

TESTS_DYNAMIC_DEPS="$CUDA_LIBRARIES $CRITERION_LIBRARIES dl m logger ${TEST_OPT_DEPS}"

if [ -n "$AMGX_FOUND" ]; then
    TESTS_DYNAMIC_DEPS="$TESTS_DYNAMIC_DEPS $AMGX_LIBRARIES"
fi

if [ -n "$SYCL_FOUND" ]; then
    TESTS_DYNAMIC_DEPS="$TESTS_DYNAMIC_DEPS mkl_sycl_blas mkl_sycl_lapack mkl_sycl_dft mkl_sycl_sparse mkl_sycl_vm mkl_sycl_rng mkl_sycl_stats mkl_sycl_data_fitting mkl_intel_ilp64 mkl_tbb_thread mkl_core sycl"
fi


##Tests
TESTS_STATIC_DEPS="monodomain ode_solver config tinyexpr config_helpers alg graph utils sds"
COMPILE_EXECUTABLE "TestSolvers" "test_solvers.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $AMGX_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH" 


TESTS_STATIC_DEPS="config alg graph utils sds"
COMPILE_EXECUTABLE "TestMesh" "test_mesh.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH" 


TESTS_STATIC_DEPS="monodomain ode_solver ini_parser config tinyexpr config_helpers alg graph utils sds"

COMPILE_OBJECT "common.c" "common.o"
COMPILE_EXECUTABLE "TestSimulations" "test_simulations.c  common.o" "common.h" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH" 


TESTS_STATIC_DEPS="alg utils sds"
COMPILE_EXECUTABLE "TestLibs" "test_libs.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH" 


##Profilers
TESTS_STATIC_DEPS="monodomain ode_solver config tinyexpr vtk_utils config_helpers alg graph utils sds miniz yxml"
COMPILE_EXECUTABLE "MeshProfiler" "profile_mesh_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

COMPILE_EXECUTABLE "VtuProfiler"  "profile_mesh_vtu_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

COMPILE_EXECUTABLE "AlgProfiler"  "profile_mesh_alg_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

COMPILE_EXECUTABLE "EnProfiler"  "profile_mesh_en_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

COMPILE_EXECUTABLE "TxtProfiler"  "profile_mesh_txt_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

COMPILE_EXECUTABLE "BinProfiler"  "profile_mesh_bin_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

COMPILE_EXECUTABLE "Custom_meshProfiler" "profile_custom_mesh_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

COMPILE_EXECUTABLE "SolversProfiler" "profile_linear_system_solvers.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

if [ -n "$CUDA_FOUND" ]; then
    COMPILE_EXECUTABLE "GpusolversProfiler" "profile_linear_system_solvers_gpu.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

fi

TESTS_STATIC_DEPS="monodomain ode_solver ini_parser config tinyexpr config_helpers alg graph utils sds"
COMPILE_EXECUTABLE "SimulationProfiler" "profile_simulation.c common.o" "common.h" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" 

