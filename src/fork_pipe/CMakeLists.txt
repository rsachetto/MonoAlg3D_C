
IF(MSVC)
    set(SOURCE_FILES_FORK  fork_pipe_windows.c fork_pipe_windows.h fork_pipe.h)
ELSE()
    set(SOURCE_FILES_FORK  fork_pipe_linux.c fork_pipe_linux.h fork_pipe.h)
ENDIF(MSVC)

ADD_LIBRARY(fork_pipe STATIC ${SOURCE_FILES_FORK})

if(NOT MSVC)
    set_property(TARGET fork_pipe PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
