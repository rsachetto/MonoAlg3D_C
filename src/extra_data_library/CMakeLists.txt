SET(LIB_TYPE SHARED) #set the lib type
SET(COMMON_SRC ../libraries_common/config_helpers.c ../libraries_common/config_helpers.c)


############## DEFAULT STIM LIBRARY 2006 ##############################
ADD_LIBRARY(default_extra_data ${LIB_TYPE} ${COMMON_SRC} extra_data.c)
#########################################################

