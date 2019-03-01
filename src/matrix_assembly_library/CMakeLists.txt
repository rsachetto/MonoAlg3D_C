SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT MATRIX ASSEMBLY LIBRARY ##############################
ADD_LIBRARY(default_matrix_assembly ${LIB_TYPE} ${COMMON_SRC} matrix_assembly.c)
TARGET_LINK_LIBRARIES(default_matrix_assembly  utils alg)
#########################################################

