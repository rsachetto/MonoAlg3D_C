SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT STIM LIBRARY ##############################
ADD_LIBRARY(default_stimuli ${LIB_TYPE} stimuli.c ${COMMON_SRC})
TARGET_LINK_LIBRARIES(default_stimuli utils alg)
#########################################################

