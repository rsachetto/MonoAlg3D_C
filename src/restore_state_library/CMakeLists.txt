SET(LIB_TYPE MODULE) #set the lib type

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(default_restore_state ${LIB_TYPE} restore_state.c)
    set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
    TARGET_LINK_LIBRARIES(default_restore_state PUBLIC  string alg)
 #   //IF(MSVC)

#ENDIF()
ELSE()
    ADD_LIBRARY(default_restore_state ${LIB_TYPE} restore_state.c)
    TARGET_LINK_LIBRARIES(default_restore_state  string alg)
ENDIF(CUDA_FOUND)

