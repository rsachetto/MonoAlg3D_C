UTILS_SOURCE_FILES="search.c stop_watch.c sort.c file_utils.c batch_utils.c heap.c"
UTILS_HEADER_FILES="utils.h stop_watch.h file_utils.h batch_utils.c heap.h"

COMPILE_STATIC_LIB "utils" "$UTILS_SOURCE_FILES" "$UTILS_HEADER_FILES"
