list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

FIND_PACKAGE(Criterion)

if(CRITERION_FOUND)
    set(SOURCE_FILES_TEST teste_cg.c)
    add_executable(TestSuite ${SOURCE_FILES_TEST})
    ADD_DEPENDENCIES(TestSuite CudaLib utils alg solvers)
    target_link_libraries(TestSuite alg solvers utils CudaLib ${CMAKE_THREAD_LIBS_INIT} m criterion)
endif()