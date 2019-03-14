SET(LIB_TYPE STATIC) #set the lib type

############## DEFAULT SAVE MESH LIBRARY ##############################

FIND_PACKAGE(ZLIB)

ADD_LIBRARY(vtk_utils ${LIB_TYPE} data_utils.c data_utils.h vtk_polydata_grid.c vtk_polydata_grid.h vtk_unstructured_grid.c vtk_unstructured_grid.h)

IF(ZLIB_FOUND)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I \"${ZLIB_INCLUDE_DIRS}\" -DCOMPILE_ZLIB")
    TARGET_LINK_LIBRARIES(vtk_utils string ${ZLIB_LIBRARIES})
ELSE()
    TARGET_LINK_LIBRARIES(vtk_utils string)
ENDIF(ZLIB_FOUND)

if(NOT MSVC)
    set_property(TARGET vtk_utils PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

#######################################################################

