SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT SAVE MESH LIBRARY ##############################

FIND_PACKAGE(ZLIB)

ADD_LIBRARY(default_save_mesh ${LIB_TYPE} ${COMMON_SRC} save_mesh.c)

IF(ZLIB_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I \"${ZLIB_INCLUDE_DIRS}\" -DCOMPILE_ZLIB")
    TARGET_LINK_LIBRARIES(default_save_mesh vtk_utils utils alg string ${ZLIB_LIBRARIES})
ELSE()
    TARGET_LINK_LIBRARIES(default_save_mesh vtk_utils utils alg string)
ENDIF(ZLIB_FOUND)

#######################################################################


