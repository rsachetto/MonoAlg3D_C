set(SOURCE_FILES_UTILS
        opts.h utils.h search.c
        point_hash.c point_hash.h
        output_utils.h output_utils.c
        ../solvers/constants.h
        stop_watch.h stop_watch.c
        string/sds.c
        string/sds.h
        string/sdsalloc.h
        sort.c
        vector/int_vector.c vector/int_vector.h
        vector/uint32_vector.c vector/uint32_vector.h)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})
