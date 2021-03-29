
if [ -n "$CUDA_FOUND" ]; then
    TEST_OPT_DEPS=gpu_utils
fi
	
	TESTS_DYNAMIC_DEPS="$CUDA_LIBRARIES $CRITERION_LIBRARIES dl m logger"

##Tests
TESTS_STATIC_DEPS="monodomain ode_solver config tinyexpr ${TEST_OPT_DEPS} config_helpers alg graph utils sds"
COMPILE_EXECUTABLE "TestSolvers" "test_solvers.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH"

TESTS_STATIC_DEPS="config alg graph utils sds"
COMPILE_EXECUTABLE "TestMesh" "test_mesh.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH"

TESTS_STATIC_DEPS="monodomain ode_solver ini_parser config tinyexpr ${TEST_OPT_DEPS} config_helpers alg graph utils sds"
COMPILE_EXECUTABLE "TestSimulations" "common.c test_simulations.c" "common.h" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH"

TESTS_STATIC_DEPS="alg utils sds"
COMPILE_EXECUTABLE "TestLibs" "test_libs.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH"

##Profilers
TESTS_STATIC_DEPS="monodomain ode_solver config tinyexpr ${TEST_OPT_DEPS} config_helpers alg graph utils sds"
COMPILE_EXECUTABLE "MeshProfiler" "profile_mesh_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY"
COMPILE_EXECUTABLE "Custom_meshProfiler" "profile_custom_mesh_load.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY"
COMPILE_EXECUTABLE "SolversProfiler" "profile_linear_system_solvers.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY"

TESTS_STATIC_DEPS="monodomain ode_solver ini_parser config tinyexpr ${TEST_OPT_DEPS} config_helpers alg graph utils sds"
COMPILE_EXECUTABLE "SimulationProfiler" "common.c profile_simulation.c" "common.h" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY"
