set(headers dlfcn.h)
set(sources dlfcn.c)
add_definitions(-DSHARED)

add_library(dl ${sources})
target_link_libraries(dl psapi)