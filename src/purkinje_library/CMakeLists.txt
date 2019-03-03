SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT PURKINJE LIBRARY ##############################
ADD_LIBRARY(default_purkinje ${LIB_TYPE} ${COMMON_SRC} purkinje.c purkinje_helpers.c purkinje_helpers.h)
TARGET_LINK_LIBRARIES(default_purkinje utils alg  graph)
######################################################################

