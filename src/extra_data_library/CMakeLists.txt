SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT STIM LIBRARY 2006 ##############################
ADD_LIBRARY(default_extra_data ${LIB_TYPE} ${COMMON_SRC} extra_data.c)
TARGET_LINK_LIBRARIES(default_extra_data hashes)
#########################################################

