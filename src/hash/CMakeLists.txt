set(SOURCE_FILES_HASH hash_common.h
        point_hash.c
        point_hash.h
        string_hash.c
        string_hash.h
        point_voidp_hash.c
        point_voidp_hash.h
        stim_config_hash.c
        stim_config_hash.h
        hash_common.h
        hash_common.c
       )

ADD_LIBRARY(hashes STATIC ${SOURCE_FILES_HASH})

if(NOT MSVC)
    set_property(TARGET hashes PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()