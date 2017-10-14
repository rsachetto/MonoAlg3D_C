set(SOURCE_FILES_STRING sds.c
        sds.h
        sdsalloc.h)

ADD_LIBRARY(string STATIC ${SOURCE_FILES_STRING})
