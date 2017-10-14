set(SOURCE_FILES_UTILS
        utils.h
        ini_parser/ini_file_sections.h
        ini_parser/ini.c
        ini_parser/ini.h
        search.c
        hash/hash_common.h
        hash/point_hash.c
        hash/point_hash.h
        hash/string_hash.c
        hash/string_hash.h
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
