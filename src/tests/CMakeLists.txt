list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

FIND_PACKAGE(Criterion)

if(CRITERION_FOUND)
    set(SOURCE_FILES_TEST tests.c)
    add_executable(TestSuite ${SOURCE_FILES_TEST})
    ADD_DEPENDENCIES(TestSuite utils alg solvers string config gpu_utils)
    target_link_libraries(TestSuite alg solvers utils string ini_parser config gpu_utils ${CMAKE_THREAD_LIBS_INIT} m criterion graph ${CUDA_LIBRARIES} )
endif()

set(SOURCE_FILES_TEST valgring_toy_test.c)
add_executable(ToyTest ${SOURCE_FILES_TEST})
ADD_DEPENDENCIES(ToyTest utils alg solvers string config gpu_utils)
target_link_libraries(ToyTest alg solvers utils string ini_parser config gpu_utils ${CMAKE_THREAD_LIBS_INIT} m )
