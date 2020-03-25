PRINT_INFO "COMPILING STATIC LIB raylib"

if [ "$BUILD_TYPE" == "debug" ]; then
  MAKE_FLAGS="RAYLIB_BUILD_MODE=DEBUG"
fi

make ${MAKE_FLAGS}

CUR_DIR=$(pwd)
COMPILED_STATIC_LIBS["raylib"]="$CUR_DIR/libraylib.a"