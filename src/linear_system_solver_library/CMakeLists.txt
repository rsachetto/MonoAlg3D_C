SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT MATRIX ASSEMBLY LIBRARY ##############################
ADD_LIBRARY(default_linear_system_solver ${LIB_TYPE} ${COMMON_SRC} linear_system_solver.c)
TARGET_LINK_LIBRARIES(default_linear_system_solver hashes utils)
#########################################################

