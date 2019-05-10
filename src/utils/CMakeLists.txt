set(SOURCE_FILES_UTILS
        utils.h
        search.c
        stop_watch.h
        stop_watch.c
        sort.c
        file_utils.c file_utils.h)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})

set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)
