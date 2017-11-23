SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT STIM LIBRARY ##############################
ADD_LIBRARY(default_stimuli ${LIB_TYPE} ${COMMON_SRC} stimuli.c)
TARGET_LINK_LIBRARIES(default_stimuli hashes utils)
#########################################################

