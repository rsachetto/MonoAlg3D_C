SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT SAVE MESH LIBRARY ##############################
ADD_LIBRARY(default_restore_state ${LIB_TYPE} ${COMMON_SRC} restore_state.c)
TARGET_LINK_LIBRARIES(default_restore_state alg)
#########################################################

