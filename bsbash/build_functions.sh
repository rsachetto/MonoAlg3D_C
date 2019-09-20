#!/bin/bash
set -e 

WARN='\033[0;31m'
INFO='\033[0;34m'
NC='\033[0m' # No Color

#MAIN BUILD VARIABLES
C_COMPILER=gcc
CXX_COMPILER=g++
AR=ar
RANLIB=ranlib
C_FLAGS=
ROOT_DIR=$(pwd)
BUILD_TYPE="release"

RUNTIME_OUTPUT_DIRECTORY=$ROOT_DIR
LIBRARY_OUTPUT_DIRECTORY=$ROOT_DIR

GLOBAL_FORCE_COMPILATION=""
QUIET=''

DEFAULT_BUILD_DIR="build_"

declare -A COMPILED_STATIC_LIBS
declare -A COMPILED_SHARED_LIBS

PRINT_USAGE () { 
	echo "Usage $0 [flags] [build_type]" >&2;
	echo "Valid flags:" >&2;
	echo "-f  - force recompilation" >&2;
	echo "-q  - quiet compilation. Only errors and warnings will be outputed" >&2;
	exit 1
}

GET_BUILD_OPTIONS () {
	
	if [ $# -gt 2 ]; then
		PRINT_USAGE "$@"
	fi

	while [ $# -gt 0 ] && [ "$1" != "--" ]; do
		while getopts "fq" opt; do
			case $opt in
			f) GLOBAL_FORCE_COMPILATION='y' ;;
			q) QUIET='y' ;;
			\?) PRINT_USAGE "$@";;
		esac
	done
		shift $((OPTIND-1))

		while [ $# -gt 0 ] && ! [[ "$1" =~ ^- ]]; do
			BUILD_TYPE=("${mandatory[@]}" "$1")
			shift
		done
	done

	if [ "$BUILD_TYPE" == "release" ]; then 
		C_FLAGS="$C_FLAGS -O3"
	elif [ "$BUILD_TYPE" == "debug" ]; then
		C_FLAGS="$C_FLAGS -g -DDEBUG_INFO"
	elif [ "$BUILD_TYPE" == "clean" ]; then
		CLEAN_PROJECT
		exit 0
	fi

}

CLEAN_PROJECT () {

	local DIR_NAME

	for T in debug release; do
		DIR_NAME="${DEFAULT_BUILD_DIR}${T}"

		if [ -d "${DIR_NAME}" ]; then
			find "${DIR_NAME}" -name "*.o" -exec rm -rf {} \;
			find "${DIR_NAME}" -name "*.a" -exec rm -rf {} \;
			find "${DIR_NAME}" -name ".*last_compiled_time_bbash" -exec rm -rf {} \;
		fi	
	done

}

PRINT_INFO () {
	if [ -z "$QUIET" ]; then
  		printf "[INFO] ${INFO}%s${NC}\n" "$1"	
	fi
}

PRINT_WARN () {
  printf "[WARN] ${WARN}%s${NC}\n" "$1" >&2
}

RECOMPILE_OR_NOT () {

  diff_time=1

  if [ -f "$1" ] && [ -f "$2" ]; then
    OBJ_TIME=$(date +%s -d"$( ls -l --time-style=full-iso $1 | awk ' { print $(NF-3)" "$(NF-2); }')")
    SOURCE_TIME=$(date +%s -d"$( ls -l --time-style=full-iso $2 | awk ' { print $(NF-3)" "$(NF-2); }')")
    (( diff_time =  SOURCE_TIME - OBJ_TIME ))
  fi

  echo $diff_time

}

ECHO_AND_EXEC_COMMAND() {
	if [ -z "$QUIET" ]; then
		echo -e "$1"
	fi
	$1
}

COMPILE_EXECUTABLE () {
  	local EXECUTABLE_NAME=$1
	local SOURCES=$2
	local HEADERS=$3
	local STATIC_DEPS_LIST=$4
	local DYNAMIC_DEPS_LIST=$5
	local EXTRA_LIB_PATH_LIST=$6
	local EXTRA_C_FLAGS=$7

	local STATIC_DEPS=()

	for dep in $STATIC_DEPS_LIST; do
	  STATIC_DEPS+=("${COMPILED_STATIC_LIBS[$dep]}")
	done

  	local DYNAMIC_DEPS=()
	for dep in $DYNAMIC_DEPS_LIST; do
	  DYNAMIC_DEPS+=("-l${dep}")
	done

  	local EXTRA_LIBRARY_PATH=()
	for libpath in $EXTRA_LIB_PATH_LIST; do
	  EXTRA_LIBRARY_PATH+=("-L${libpath}")
	done

  	local BUILD_DIR=$ROOT_DIR/$RUNTIME_OUTPUT_DIRECTORY

	if [ ! -d "$BUILD_DIR" ]; then
		mkdir -p "$BUILD_DIR"
	fi

  	PRINT_INFO "COMPILING AND LINKING EXECUTABLE $EXECUTABLE_NAME"

	local FORCE_COMPILATION=$GLOBAL_FORCE_COMPILATION

	local TIME_FILE="$BUILD_DIR/.${EXECUTABLE_NAME}_last_compiled_time_bbash"

	if [ -z "$FORCE_COMPILATION" ]; then
		#CHECK IF BUILD SCRIPT CHANGED SINCE THE LAST COMPILATION
		if [ "$(RECOMPILE_OR_NOT "$TIME_FILE" "$0")" -gt "0" ]; then
			FORCE_COMPILATION='y'
		fi
	fi
	
	if [ -z "$FORCE_COMPILATION" ]; then
		for h in $HEADERS; do
			if [ -f "$h" ]; then
			if [ "$(RECOMPILE_OR_NOT "$BUILD_DIR/${EXECUTABLE_NAME}" "$h")" -gt "0" ]; then
				FORCE_COMPILATION='y'
				break
			fi
			else
				PRINT_WARN "$h file does not exist! Check your build scripts!"
			fi
		done
	fi

	if [ -z "$FORCE_COMPILATION" ]; then
		for s in $SOURCES; do
			if [ -f "$s" ]; then
				if [ "$(RECOMPILE_OR_NOT "$BUILD_DIR/${EXECUTABLE_NAME}" "$s")" -gt "0" ]; then
					FORCE_COMPILATION='y'
					break
				fi
			fi
		done
	fi

	if [ -z "$FORCE_COMPILATION" ]; then
		for d in $SOURCES; do
			if [ -f "$s" ]; then
				if [ "$(RECOMPILE_OR_NOT "$BUILD_DIR/${EXECUTABLE_NAME}" "{STATIC_DEPS[$d]}")" -gt "0" ]; then
					FORCE_COMPILATION='y'
					break
				fi
			fi
		done
	fi

	if [ -n "$FORCE_COMPILATION" ]; then
		local MY_C_FLAGS="$C_FLAGS $EXTRA_C_FLAGS"
		ECHO_AND_EXEC_COMMAND "gcc $MY_C_FLAGS $SOURCES ${STATIC_DEPS[*]} -o $BUILD_DIR/${EXECUTABLE_NAME} ${EXTRA_LIBRARY_PATH[*]} ${DYNAMIC_DEPS[*]}"
		touch $TIME_FILE
	fi

}

COMPILE_OBJECT () {

	local SRC_FILE=$1
	local OBJ_FILE=$2

	local MY_C_FLAGS="$C_FLAGS $3"
	local FORCE_COMPILATION=$4
	local IS_CUDA=$5

	local COMPILER=$C_COMPILER

	PRINT_INFO "COMPILING OBJECT $OBJ_FILE"

	ANY_COMPILED=""

	if [ -n "$IS_CUDA" ] && [[ $SRC_FILE == *.cu ]]; then
		COMPILER=$NVCC
	elif [[ $SRC_FILE == *.cpp ]] || [[ $SRC_FILE == *.cxx ]]; then
		COMPILER=$CXX_COMPILER
	fi

	if [ -n "$FORCE_COMPILATION" ] || [ "$(RECOMPILE_OR_NOT "$OBJ_FILE" "$SRC_FILE")" -gt "0" ]; then
		if [[ "$COMPILER" == "$NVCC" ]]; then
			local X_COMPILER_FLAGS=()
			for xflag in $MY_C_FLAGS; do
			if [ "$xflag" != "-std=gnu99" ]; then
				X_COMPILER_FLAGS+=("\\\"${xflag}\\\",")
			fi
			done

			local S_X_COMPILER_FLAGS
			S_X_COMPILER_FLAGS=$( printf "%s" "${X_COMPILER_FLAGS[@]}" )

			ECHO_AND_EXEC_COMMAND "$COMPILER $SRC_FILE -c  -o $OBJ_FILE -ccbin $C_COMPILER -m64 -Xcompiler ${S_X_COMPILER_FLAGS} -DNVCC -I$CUDA_INCLUDE_PATH"

		else 
			ECHO_AND_EXEC_COMMAND "$COMPILER $MY_C_FLAGS -c $SRC_FILE -o $OBJ_FILE"
		fi

		ANY_COMPILED="y"

	fi

}

COMPILE_STATIC_LIB () {
 
	local LIB_NAME=lib$1
	local SOURCES=$2
	local HEADERS=$3
	local EXTRA_C_FLAGS=$4

	local BUILD_DIR=$ROOT_DIR/${DEFAULT_BUILD_DIR}${BUILD_TYPE}/$LIB_NAME
	local LIB_PATH=$BUILD_DIR/$LIB_NAME.a

	if [ ! -d "$BUILD_DIR/objs" ]; then
		mkdir -p "$BUILD_DIR/objs"
	fi

	local OBJECTS=()

	PRINT_INFO "COMPILING STATIC LIB $LIB_NAME"

	local ANY_COMPILED_LOCAL=''

	local FORCE_COMPILATION=$GLOBAL_FORCE_COMPILATION

	local TIME_FILE="$BUILD_DIR/.${LIB_NAME}_last_compiled_time_bbash"

	
	if [ -z "$FORCE_COMPILATION" ]; then

		#CHECK IF BUILD SCRIPT CHANGED SINCE THE LAST COMPILATION
		if [ "$(RECOMPILE_OR_NOT "$TIME_FILE" "$0")" -gt "0" ]; then
			FORCE_COMPILATION='y'
		fi


		for h in $HEADERS; do
			if [ -f "$h" ]; then
			if [ "$(RECOMPILE_OR_NOT "$LIB_PATH" "$h")" -gt "0" ]; then
				FORCE_COMPILATION='y'
				break
			fi
			else
			PRINT_WARN "$h file does not exist! Check your build scripts!"
			fi
		done
	
	fi

	for s in $SOURCES; do
		local OBJ_FILE
		OBJ_FILE=$(basename -- "$s")
		OBJ_FILE=$BUILD_DIR/objs/${OBJ_FILE}.o
		OBJECTS+=("$OBJ_FILE")

		COMPILE_OBJECT "$s" "$OBJ_FILE" "$EXTRA_C_FLAGS -fPIC" "$FORCE_COMPILATION"

		if [ -z "$ANY_COMPILED_LOCAL" ]; then
		ANY_COMPILED_LOCAL=$ANY_COMPILED
		fi

	done

  	PRINT_INFO "CREATING STATIC LIB $LIB_NAME"

	if [ -n "$ANY_COMPILED_LOCAL" ]; then
		ECHO_AND_EXEC_COMMAND "$AR rcs $LIB_PATH ${OBJECTS[*]}"
		ECHO_AND_EXEC_COMMAND "$RANLIB $LIB_PATH"
		touch $TIME_FILE
	fi

	COMPILED_STATIC_LIBS[$1]=$LIB_PATH

}

COMPILE_SHARED_LIB () {

	local LIB_NAME=lib$1
	local SOURCES=$2
	local HEADERS=$3
	local STATIC_DEPS_LIST=$4
	local DYNAMIC_DEPS_LIST=$5
	local EXTRA_LIB_PATH_LIST=$6
	local EXTRA_C_FLAGS=$7
	local IS_CUDA=$8

	local STATIC_DEPS=()

	for dep in $STATIC_DEPS_LIST; do
	  STATIC_DEPS+=("${COMPILED_STATIC_LIBS[$dep]}")
	done

  	local DYNAMIC_DEPS=()
	for dep in $DYNAMIC_DEPS_LIST; do
	  DYNAMIC_DEPS+=("-l${dep}")
	done

  	local EXTRA_LIBRARY_PATH=()
	for libpath in $EXTRA_LIB_PATH_LIST; do
	  if [ -n "${libpath}" ]; then
	    EXTRA_LIBRARY_PATH+=("-L${libpath}")
	  fi
	done

	local BUILD_DIR=$ROOT_DIR/${DEFAULT_BUILD_DIR}${BUILD_TYPE}/$LIB_NAME

  	if [ ! -d "$BUILD_DIR/objs" ]; then
		mkdir -p "$BUILD_DIR/objs"
	fi

	local OBJECTS=()

	local LIB_PATH=$BUILD_DIR/$LIB_NAME.so

	PRINT_INFO "COMPILING SHARED LIB $LIB_NAME"

	local ANY_COMPILED_LOCAL=''

   	local FORCE_COMPILATION=$GLOBAL_FORCE_COMPILATION
	local TIME_FILE="$BUILD_DIR/.${LIB_NAME}_last_compiled_time_bbash"

	if [ -z $FORCE_COMPILATION ]; then
		#CHECK IF BUILD SCRIPT CHANGED SINCE THE LAST COMPILATION
		if [ "$(RECOMPILE_OR_NOT "$TIME_FILE" "$0")" -gt "0" ]; then
			FORCE_COMPILATION='y'
		fi

		for h in $HEADERS; do
			if [ -f "$h" ]; then
			if [ "$(RECOMPILE_OR_NOT "$LIB_PATH" "$h")" -gt "0" ]; then
				FORCE_COMPILATION='y'
				break
			fi
			else
				PRINT_WARN "$h file does not exist! Check your build scripts!"
			fi
		done
	fi

	for s in $SOURCES; do
		local OBJ_FILE
		OBJ_FILE=$(basename -- "$s")
		OBJ_FILE=$BUILD_DIR/objs/${OBJ_FILE}.o
		OBJECTS+=("$OBJ_FILE")

    	COMPILE_OBJECT "$s" "$OBJ_FILE" "${EXTRA_C_FLAGS} -fPIC" "$FORCE_COMPILATION" "$IS_CUDA"

		if [ -z "$ANY_COMPILED_LOCAL" ]; then
			ANY_COMPILED_LOCAL=$ANY_COMPILED
		fi

	done

  	PRINT_INFO "LINKING SHARED LIB $LIB_NAME"

	if [ -n "$IS_CUDA" ]; then
		LINKER=$CXX_COMPILER
	else
		LINKER=$C_COMPILER
	fi

	if [ -n "$ANY_COMPILED_LOCAL" ]; then
		ALL_FLAGS="-fPIC $C_FLAGS -shared -o $LIB_PATH ${OBJECTS[*]} ${STATIC_DEPS[*]} ${EXTRA_LIBRARY_PATH[*]} ${DYNAMIC_DEPS[*]}"
		if [ -n "$IS_CUDA" ]; then
			ECHO_AND_EXEC_COMMAND "$LINKER $ALL_FLAGS"
		else
			ECHO_AND_EXEC_COMMAND "$LINKER $ALL_FLAGS"
		fi
		
		ECHO_AND_EXEC_COMMAND "cp $LIB_PATH $ROOT_DIR/$LIBRARY_OUTPUT_DIRECTORY"
		touch $TIME_FILE

	fi

  	COMPILED_SHARED_LIBS[$1]=$LIB_PATH

}

ADD_SUBDIRECTORY() {
	local PREVIOUS_DIR
	PREVIOUS_DIR=$(pwd)
	cd "$1" || exit
	source build.sh
	cd "$PREVIOUS_DIR" || exit
}