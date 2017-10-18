set(SOURCE_FILES_UTILS
        utils.h
        search.c
        stop_watch.h
        stop_watch.c
        sort.c erros_helpers.c erros_helpers.h)

ADD_LIBRARY(utils SHARED ${SOURCE_FILES_UTILS})
