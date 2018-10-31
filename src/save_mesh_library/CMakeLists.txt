SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT SAVE MESH LIBRARY ##############################
ADD_LIBRARY(default_save_mesh ${LIB_TYPE} ${COMMON_SRC} save_mesh.c vtk_unstructured_grid.c vtk_unstructured_grid.h)
TARGET_LINK_LIBRARIES(default_save_mesh hashes utils alg string)
#########################################################

