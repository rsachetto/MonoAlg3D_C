SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT MODIFY DOMAIN LIBRARY ##############################
ADD_LIBRARY(default_modify_domain ${LIB_TYPE} modify_domain.c ${COMMON_SRC})
TARGET_LINK_LIBRARIES(default_modify_domain utils alg)
#########################################################

