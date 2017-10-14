set(SOURCE_FILES_UTILS
        utils.h
        search.c
        stop_watch.h
        stop_watch.c
        sort.c)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})
