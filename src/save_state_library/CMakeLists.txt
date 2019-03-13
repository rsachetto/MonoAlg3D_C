SET(LIB_TYPE MODULE) #set the lib type
if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(default_save_state ${LIB_TYPE}  save_state.c)
    set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
    TARGET_LINK_LIBRARIES(default_save_state PUBLIC  string alg)
ELSE()
    ADD_LIBRARY(default_save_state ${LIB_TYPE} save_state.c)
    TARGET_LINK_LIBRARIES(default_save_state string alg)
ENDIF(CUDA_FOUND)




