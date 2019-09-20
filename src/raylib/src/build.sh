PRINT_INFO "COMPILING STATIC LIB raylib"

make
CUR_DIR=$(pwd)
COMPILED_STATIC_LIBS["raylib"]="$CUR_DIR/libraylib.a"