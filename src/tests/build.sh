
if [ -n "$CUDA_FOUND" ]; then
    TEST_OPT_DEPS=gpu_utils
fi

if [ -n "$COMPILE_GUI" ]; then
	TESTS_STATIC_DEPS="monodomain ode_solver ini_parser config tinyexpr ${TEST_OPT_DEPS} utils config_helpers  alg graph utils sds gui vtk_utils miniz yxml raylib tinyfiledialogs"
	TESTS_DYNAMIC_DEPS="$CUDA_LIBRARIES $CRITERION_LIBRARIES dl m logger OpenGL GLX GLU pthread X11 rt"
else
	TESTS_STATIC_DEPS="monodomain ode_solver ini_parser config tinyexpr ${TEST_OPT_DEPS} utils config_helpers vtk_utils yxml alg graph utils sds miniz"
	TESTS_DYNAMIC_DEPS="$CUDA_LIBRARIES $CRITERION_LIBRARIES dl m logger"
fi

##Tests
COMPILE_EXECUTABLE "TestSolvers" "test_solvers.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH"
COMPILE_EXECUTABLE "TestMesh" "test_mesh.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH"
COMPILE_EXECUTABLE "TestSimulations" "test_simulations.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH"
COMPILE_EXECUTABLE "TestLibs" "test_libs.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS" "$CUDA_LIBRARY_PATH $CRITERION_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY" "$CRITERION_INCLUDE_PATH"

##Profilers
#COMPILE_OBJECT "malloc_count.c" "malloc_count.o"
#COMPILE_EXECUTABLE "MeshProfiler" "profile_mesh_load.c malloc_count.o" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY"
#COMPILE_EXECUTABLE "SolversProfiler" "profile_linear_system_solvers.c" "" "$TESTS_STATIC_DEPS" "$TESTS_DYNAMIC_DEPS gdbm" "$CUDA_LIBRARY_PATH $LIBRARY_OUTPUT_DIRECTORY"

#rm malloc_count.o
