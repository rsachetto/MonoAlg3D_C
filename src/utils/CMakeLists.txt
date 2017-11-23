set(SOURCE_FILES_UTILS
        utils.h
        search.c
        stop_watch.h
        stop_watch.c
        sort.c
        logfile_utils.c logfile_utils.h)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})

if(NOT MSVC)
    set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
