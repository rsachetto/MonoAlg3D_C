set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests_bin )

set(SOURCE_FILES_TEST teste_cg.c)

add_executable(TestSuite ${SOURCE_FILES_TEST})

ADD_DEPENDENCIES(TestSuite CudaLib utils alg solvers)

target_link_libraries(TestSuite alg solvers utils CudaLib ${CMAKE_THREAD_LIBS_INIT} m criterion)
