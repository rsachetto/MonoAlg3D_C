set(SOURCE_FILES_STRING sds.c
        sds.h
        sdsalloc.h)

ADD_LIBRARY(string STATIC ${SOURCE_FILES_STRING})

set_property(TARGET string PROPERTY POSITION_INDEPENDENT_CODE ON)
