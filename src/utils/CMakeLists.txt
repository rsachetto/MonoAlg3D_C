set(SOURCE_FILES_UTILS
        opts.h utils.h search.c
        point_hash.c point_hash.h
        output_utils.h output_utils.c
        constants.h
        stop_watch.h stop_watch.c
        string/sds.c
        string/sds.h
        string/sdsalloc.h)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})
