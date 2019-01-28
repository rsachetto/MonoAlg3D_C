SET(LIB_TYPE MODULE data_utils.c data_utils.h) #set the lib type

############## DEFAULT SAVE MESH LIBRARY ##############################

FIND_PACKAGE(ZLIB)

ADD_LIBRARY(default_save_mesh ${LIB_TYPE} ${COMMON_SRC} save_mesh.c vtk_unstructured_grid.c vtk_unstructured_grid.h vtk_polydata_grid.c vtk_polydata_grid.h data_utils.c data_utils.h)

IF(ZLIB_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I \"${ZLIB_INCLUDE_DIRS}\" -DCOMPILE_OPENGL")
    TARGET_LINK_LIBRARIES(default_save_mesh hashes utils alg string ${ZLIB_LIBRARIES})
ELSE()
    TARGET_LINK_LIBRARIES(default_save_mesh hashes utils alg string)
ENDIF(ZLIB_FOUND)

#######################################################################

