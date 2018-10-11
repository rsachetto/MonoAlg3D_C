SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT SAVE MESH LIBRARY ##############################
ADD_LIBRARY(default_save_state ${LIB_TYPE} ${COMMON_SRC} save_state.c)
TARGET_LINK_LIBRARIES(default_save_state alg)
#########################################################

