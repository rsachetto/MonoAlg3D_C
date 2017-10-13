set(SOURCE_FILES_UTILS
        utils.h
        config_parser.h
        config_parser.c
        search.c
        point_hash.c
        point_hash.h
        output_utils.h
        output_utils.c
        stop_watch.h
        stop_watch.c
        string/sds.c
        string/sds.h
        string/sdsalloc.h
        sort.c
        vector/int_vector.c vector/int_vector.h
        vector/uint32_vector.c vector/uint32_vector.h
        gpu/gpu_utils.c gpu/gpu_utils.h)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})
