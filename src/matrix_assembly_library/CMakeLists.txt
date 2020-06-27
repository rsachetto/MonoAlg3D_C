SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT MATRIX ASSEMBLY LIBRARY ####################################
ADD_LIBRARY(default_matrix_assembly ${LIB_TYPE} ${COMMON_SRC} matrix_assembly.c)
TARGET_LINK_LIBRARIES(default_matrix_assembly  utils alg config_helpers)
####################################################################################

############## DDM MATRIX ASSEMBLY LIBRARY #########################################
ADD_LIBRARY(ddm_matrix_assembly ${LIB_TYPE} ${COMMON_SRC} ddm_matrix_assembly.c)
TARGET_LINK_LIBRARIES(ddm_matrix_assembly  utils alg config_helpers)
####################################################################################

############## PURKINJE MATRIX ASSEMBLY LIBRARY ####################################
ADD_LIBRARY(purkinje_matrix_assembly ${LIB_TYPE} ${COMMON_SRC} purkinje_matrix_assembly.c)
TARGET_LINK_LIBRARIES(purkinje_matrix_assembly  utils alg config_helpers)
####################################################################################

############## PURKINJE-COUPLING MATRIX ASSEMBLY LIBRARY ####################################
ADD_LIBRARY(purkinje_coupling_matrix_assembly ${LIB_TYPE} ${COMMON_SRC} purkinje_coupling_matrix_assembly.c)
TARGET_LINK_LIBRARIES(purkinje_coupling_matrix_assembly utils alg config_helpers)
####################################################################################

