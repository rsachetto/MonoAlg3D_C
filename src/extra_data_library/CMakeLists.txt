SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT EXTRA LIBRARY ####################################
ADD_LIBRARY(default_extra_data ${LIB_TYPE} ${COMMON_SRC} extra_data.c)
TARGET_LINK_LIBRARIES(default_extra_data alg)
#########################################################################

############## MIXED CELULAR MODELS LIBRARY #############################
ADD_LIBRARY(mixed_models_extra_data ${LIB_TYPE} ${COMMON_SRC} extra_mixed_models.c)
TARGET_LINK_LIBRARIES(mixed_models_extra_data alg)
#########################################################################

