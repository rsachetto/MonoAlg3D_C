# Define main library target
add_library(yxml STATIC yxml.c)
set_property(TARGET yxml PROPERTY POSITION_INDEPENDENT_CODE ON)
