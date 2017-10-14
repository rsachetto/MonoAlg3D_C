set(SOURCE_FILES_HASH hash_common.h
        point_hash.c
        point_hash.h
        string_hash.c
        string_hash.h
       )

ADD_LIBRARY(hashes SHARED ${SOURCE_FILES_HASH})
