SET(LIB_TYPE MODULE) #set the lib type

if(CUDA_FOUND)
    set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
    ############## DEFAULT MATRIX ASSEMBLY LIBRARY ##############################
    CUDA_ADD_LIBRARY(default_linear_system_solver ${LIB_TYPE} ${COMMON_SRC} linear_system_solver.c)
    CUDA_ADD_CUBLAS_TO_TARGET(default_linear_system_solver)
    TARGET_LINK_LIBRARIES(default_linear_system_solver PUBLIC utils alg ${CUDA_cusparse_LIBRARY})
else()
    ADD_LIBRARY(default_linear_system_solver ${LIB_TYPE} ${COMMON_SRC} linear_system_solver.c)
    TARGET_LINK_LIBRARIES(default_linear_system_solver utils alg)
endif(CUDA_FOUND)
#########################################################

