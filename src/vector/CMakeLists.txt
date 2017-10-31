set(VECTOR_FILES_UTILS
        int_vector.c int_vector.h
        uint32_vector.c uint32_vector.h
        #element_vector.c
        #element_vector.h
        )

ADD_LIBRARY(vector STATIC ${VECTOR_FILES_UTILS})
