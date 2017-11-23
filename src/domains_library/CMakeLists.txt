SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT STIM LIBRARY 2006 ##############################
ADD_LIBRARY(default_domains ${LIB_TYPE} ${COMMON_SRC} domain.c domain_helpers.c domain_helpers.h)
TARGET_LINK_LIBRARIES(default_domains utils alg hashes)
#########################################################

