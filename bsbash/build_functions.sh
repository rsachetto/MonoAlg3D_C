#!/bin/bash
set -e

WARN='\033[0;31m'
INFO='\033[0;34m'
ERROR='\033[0;31m'
NC='\033[0m' # No Color

#MAIN BUILD VARIABLES
C_COMPILER="gcc"
CXX_COMPILER="g++"
AR_COMMAND="ar"
RANLIB="ranlib"
C_FLAGS=
ROOT_DIR=$(pwd)
BUILD_TYPE="release"

RUNTIME_OUTPUT_DIRECTORY=$ROOT_DIR
LIBRARY_OUTPUT_DIRECTORY=$ROOT_DIR

GLOBAL_FORCE_COMPILATION=""
QUIET=''
WAIT_ENTER=''
WRITE_COMPILE_COMMANDS='y'

DEFAULT_BUILD_DIR="build_"
COMPILE_COMMANDS_FILE="${ROOT_DIR}/compile_commands.json"

declare -A COMPILED_STATIC_LIBS
declare -A COMPILED_SHARED_LIBS

CHECK_COLOR_SUPPORT() {

    if [ -z "$TERM" ] || [ "$TERM" == "dumb" ]; then
        WARN=''
        INFO=''
        ERROR=''
        NC=''
        return
    fi

    #check if the terminal supports colored output
    command -v tput >/dev/null 2>&1 && {

        COLOR=$(tput colors)

        if [ "$COLOR" -eq -1 ]; then
            WARN=''
            INFO=''
            ERROR=''
            NC=''
        fi
    }
}

GET_LINUX_VERSION() {
	if [ -f /etc/os-release ]; then
		# freedesktop.org and systemd
		. /etc/os-release
		OS=$NAME
		VER=$VERSION_ID
	elif type lsb_release >/dev/null 2>&1; then
		# linuxbase.org
		OS=$(lsb_release -si)
		VER=$(lsb_release -sr)
	elif [ -f /etc/lsb-release ]; then
		# For some versions of Debian/Ubuntu without lsb_release command
		. /etc/lsb-release
		OS=$DISTRIB_ID
		VER=$DISTRIB_RELEASE
	elif [ -f /etc/debian_version ]; then
		# Older Debian/Ubuntu/etc.
		OS=Debian
		VER=$(cat /etc/debian_version)
	elif [ -f /etc/SuSe-release ]; then
		# Older SuSE/etc.
		...
	elif [ -f /etc/redhat-release ]; then
		# Older Red Hat, CentOS, etc.
		...
	else
		# Fall back to uname, e.g. "Linux <version>", also works for BSD, etc.
		OS=$(uname -s)
		VER=$(uname -r)
	fi

	export OS
	export VER
}

BUILD_TYPE_FILE=".build_type"

GET_BUILD_OPTIONS () {

	OPTIND=1

	while getopts "h?fqwlrd" opt; do
		case "$opt" in
			h|\?) PRINT_USAGE "$@" ;;
			f) GLOBAL_FORCE_COMPILATION='y' ;;
			l) WRITE_COMPILE_COMMANDS='y'; GLOBAL_FORCE_COMPILATION='y' ;;
			r) BUILD_TYPE='release' ;;
			d) BUILD_TYPE='debug' ;;
			q) QUIET='y' ;;
            w) WAIT_ENTER='y' ;;
		esac
	done

	shift $((OPTIND-1))

	[ "${1:-}" = "--" ] && shift

	BUILD_ARGS=( "$@" )

	if [ ${#BUILD_ARGS[@]} -eq 0 ]; then
		BUILD_ARGS+=('all')
	fi

    NEW=$BUILD_TYPE

    if [ ! -f "$BUILD_TYPE_FILE" ]; then
        LAST=$BUILD_TYPE
    else
        LAST=$(grep new .build_type | awk -F "=" '{print $2}')
    fi

    echo "last=$LAST" >  $BUILD_TYPE_FILE
    echo "new=$NEW"  >>  $BUILD_TYPE_FILE

    if [ "$NEW" != "$LAST" ]; then
       rm -f "$COMPILE_COMMANDS_FILE"
    fi

}

CLEAN_PROJECT () {

	DIR_NAME="${DEFAULT_BUILD_DIR}${1}"

	if [ -d "${DIR_NAME}" ]; then
		rm -fr "${DIR_NAME}"
	fi

}

PRINT_INFO () {
	if [ -z "$QUIET" ]; then
		printf "[INFO] ${INFO}%s${NC}\n" "$1"
	fi
}

PRINT_WARN () {
	printf "[WARN] ${WARN}%s${NC}\n" "$1" >&2
}

PRINT_ERROR () {
	printf "[ERROR] ${ERROR}%s${NC}\n" "$1" >&2
}


RECOMPILE_OR_NOT () {

	diff_time=1

	if [ -f "$1" ] && [ -f "$2" ]; then
		#shellcheck disable=SC2012
		OBJ_TIME=$(date +%s -d"$( ls -l --time-style=full-iso "$1" | awk ' { print $(NF-3)" "$(NF-2); }')")
		#shellcheck disable=SC2012
		SOURCE_TIME=$(date +%s -d"$( ls -l --time-style=full-iso "$2" | awk ' { print $(NF-3)" "$(NF-2); }')")
		(( diff_time =  SOURCE_TIME - OBJ_TIME ))
	fi

	echo "$diff_time"

}

CREATE_COMPILE_COMMANDS_FILE() {
    echo "[" > "$COMPILE_COMMANDS_FILE"
    echo "]" >> "$COMPILE_COMMANDS_FILE"
    INSERT_TO_LINE=2
}

ADD_COMPILE_COMMAND() {

	COMMAND=$1
	FILE_FULL_PATH=$2

	printf -v ESCAPED_COMMAND "%q" "$COMMAND"

	local FILE
	FILE=$(basename "$FILE_FULL_PATH")
	FILES=($FILE)

	local DIR
	DIR=$(dirname "$FILE_FULL_PATH")

    if [ ! -f "$COMPILE_COMMANDS_FILE" ]; then
	    CREATE_COMPILE_COMMANDS_FILE
	else
	    INSERT_TO_LINE=$(wc -l < "$COMPILE_COMMANDS_FILE")
	fi

    COMMAND_EXISTS=$(grep -c "$COMMAND" "$COMPILE_COMMANDS_FILE") || true #avoid exiting the script

    if [ "$COMMAND_EXISTS" == "0" ]; then
        sed -i "${INSERT_TO_LINE}i{" "$COMPILE_COMMANDS_FILE"
        INSERT_TO_LINE=$((INSERT_TO_LINE+1))

        sed -i "${INSERT_TO_LINE}i\"directory\": \"$DIR\", " "$COMPILE_COMMANDS_FILE"
        INSERT_TO_LINE=$((INSERT_TO_LINE+1))

        sed -i "${INSERT_TO_LINE}i\"command\": \"$ESCAPED_COMMAND\", " "$COMPILE_COMMANDS_FILE"
        INSERT_TO_LINE=$((INSERT_TO_LINE+1))

        sed -i "${INSERT_TO_LINE}i\"file\": \"${FILES[0]}\" " "$COMPILE_COMMANDS_FILE"
        INSERT_TO_LINE=$((INSERT_TO_LINE+1))

        sed -i "${INSERT_TO_LINE}i    }," "$COMPILE_COMMANDS_FILE"
        INSERT_TO_LINE=$((INSERT_TO_LINE+1))
	fi

}

ECHO_AND_EXEC_COMMAND() {

	if [ -z "$QUIET" ]; then
		echo -e "$1"
	fi

	#EXEC THE COMMAND
	$1

	ret_val=$?

	if [ $ret_val -ne 0 ]; then
		echo "Error issuing command $1"
		exit $ret_val
	fi

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

	local BUILD_DIR=$RUNTIME_OUTPUT_DIRECTORY

	if [ ! -d "$BUILD_DIR" ]; then
		mkdir -p "$BUILD_DIR"
	fi

	local FORCE_COMPILATION=$GLOBAL_FORCE_COMPILATION

	local TIME_FILE="$BUILD_DIR/.${EXECUTABLE_NAME}_last_compiled_time_bbash"

	#CHECK IF BUILD SCRIPT CHANGED SINCE THE LAST COMPILATION
	if [ -z "$FORCE_COMPILATION" ]; then
		if [ "$(RECOMPILE_OR_NOT "$TIME_FILE" "$0")" -gt "0" ]; then
			FORCE_COMPILATION='y'
		fi
	fi

	#CHECK IF ANY HEADER CHANGED SINCE THE LAST COMPILATION
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

	#CHECK IF ANY SOURCE CHANGED SINCE THE LAST COMPILATION
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

	#CHECK IF ANY STATIC DEPENDENCY CHANGED SINCE THE LAST COMPILATION
	if [ -z "$FORCE_COMPILATION" ]; then
		for d in ${STATIC_DEPS[*]}; do
			if [ -f "$d" ]; then
				if [ "$(RECOMPILE_OR_NOT "$BUILD_DIR/${EXECUTABLE_NAME}" "$d")" -gt "0" ]; then
					FORCE_COMPILATION='y'
					break
				fi
			fi
		done
	fi

	#CHECK IF ANY DYNAMIC DEPENDENCY CHANGED SINCE THE LAST COMPILATION
	if [ -z "$FORCE_COMPILATION" ]; then
		for d in ${DYNAMIC_DEPS[*]}; do
			if [ -f "$d" ]; then
				if [ "$(RECOMPILE_OR_NOT "$BUILD_DIR/${EXECUTABLE_NAME}" "$d")" -gt "0" ]; then
					FORCE_COMPILATION='y'
					break
				fi
			fi
		done
	fi

	if [ -n "$FORCE_COMPILATION" ]; then

		PRINT_INFO "COMPILING AND LINKING EXECUTABLE $EXECUTABLE_NAME"

		local MY_C_FLAGS="$C_FLAGS $EXTRA_C_FLAGS"

		local COMPILER_COMMAND="gcc $MY_C_FLAGS $SOURCES ${STATIC_DEPS[*]}  -o $BUILD_DIR/${EXECUTABLE_NAME} ${EXTRA_LIBRARY_PATH[*]} ${DYNAMIC_DEPS[*]} -Wl,-rpath=$LIBRARY_OUTPUT_DIRECTORY"

		ECHO_AND_EXEC_COMMAND "${COMPILER_COMMAND}"
		touch "$TIME_FILE"

		if [ -n "$WRITE_COMPILE_COMMANDS" ] ; then
			ADD_COMPILE_COMMAND "$COMPILER_COMMAND"  "$PWD/${SOURCES}"
		fi

	fi
}

COMPILE_OBJECT () {

	local SRC_FILE=$1
	local OBJ_FILE=$2

	local MY_C_FLAGS="$C_FLAGS $3"
	local FORCE_COMPILATION=$4
	local IS_CUDA=$5

	local COMPILER=$C_COMPILER

	ANY_COMPILED=""

	if [ -n "$IS_CUDA" ] && [[ $SRC_FILE == *.cu ]]; then
		COMPILER=$NVCC
	elif [[ $SRC_FILE == *.cpp ]] || [[ $SRC_FILE == *.cxx ]]; then
		COMPILER=$CXX_COMPILER
	fi

	COMPILER_COMMAND=''

	if [ -n "$FORCE_COMPILATION" ] || [ "$(RECOMPILE_OR_NOT "$OBJ_FILE" "$SRC_FILE")" -gt "0" ]; then
		if [[ "$COMPILER" == "$NVCC" ]]; then
			local X_COMPILER_FLAGS=()
			for xflag in $MY_C_FLAGS; do
				if [ "$xflag" != "-std=gnu99" ] && [ "$xflag" != "-Werror=implicit-function-declaration" ]; then
					X_COMPILER_FLAGS+=("\\\"${xflag}\\\",")
				fi
			done

			local S_X_COMPILER_FLAGS
			S_X_COMPILER_FLAGS=$( printf "%s" "${X_COMPILER_FLAGS[@]}" )
			COMPILER_COMMAND="$COMPILER $SRC_FILE -c  -o $OBJ_FILE -ccbin $C_COMPILER -m64 -Xcompiler ${S_X_COMPILER_FLAGS} -DNVCC -I$CUDA_INCLUDE_PATH"

		else
			COMPILER_COMMAND="$COMPILER $MY_C_FLAGS -c $SRC_FILE -o $OBJ_FILE"
			if [ -n "$WRITE_COMPILE_COMMANDS" ] ; then
				ADD_COMPILE_COMMAND "$COMPILER_COMMAND" "$SRC_FILE"
			fi

		fi

		PRINT_INFO "COMPILING OBJECT $OBJ_FILE"
		ECHO_AND_EXEC_COMMAND "$COMPILER_COMMAND"

		ANY_COMPILED="y"

	fi



}

CHECK_HEADERS() {
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

	local ANY_COMPILED_LOCAL=''

	local FORCE_COMPILATION=$GLOBAL_FORCE_COMPILATION

	local TIME_FILE="$BUILD_DIR/.${LIB_NAME}_last_compiled_time_bbash"

	CHECK_HEADERS

	for s in $SOURCES; do
		local OBJ_FILE
		OBJ_FILE=$(basename -- "$s")
		OBJ_FILE=$BUILD_DIR/objs/${OBJ_FILE}.o
		OBJECTS+=("$OBJ_FILE")

		COMPILE_OBJECT "${PWD}/$s" "$OBJ_FILE" "$EXTRA_C_FLAGS -fPIC" "$FORCE_COMPILATION"

		if [ -z "$ANY_COMPILED_LOCAL" ]; then
			ANY_COMPILED_LOCAL=$ANY_COMPILED
		fi
	done

	if [ -n "$ANY_COMPILED_LOCAL" ]; then
		PRINT_INFO "CREATING STATIC LIB $LIB_NAME"
		ECHO_AND_EXEC_COMMAND "$AR_COMMAND rcs $LIB_PATH ${OBJECTS[*]}"
		ECHO_AND_EXEC_COMMAND "$RANLIB $LIB_PATH"
		touch "$TIME_FILE"
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

	local ANY_COMPILED_LOCAL=''

	local FORCE_COMPILATION=$GLOBAL_FORCE_COMPILATION
	local TIME_FILE="$BUILD_DIR/.${LIB_NAME}_last_compiled_time_bbash"

	if [ ! -d "$LIBRARY_OUTPUT_DIRECTORY" ]; then
		mkdir -p "$LIBRARY_OUTPUT_DIRECTORY"
		FORCE_COMPILATION='y'
	fi

	CHECK_HEADERS

	for s in $SOURCES; do
		local OBJ_FILE
		OBJ_FILE=$(basename -- "$s")
		OBJ_FILE=$BUILD_DIR/objs/${OBJ_FILE}.o
		OBJECTS+=("$OBJ_FILE")

		COMPILE_OBJECT "${PWD}/$s" "$OBJ_FILE" "${EXTRA_C_FLAGS} -fPIC" "$FORCE_COMPILATION" "$IS_CUDA"

		if [ -z "$ANY_COMPILED_LOCAL" ]; then
			ANY_COMPILED_LOCAL=$ANY_COMPILED
		fi
	done

	if [ -n "$IS_CUDA" ]; then
		LINKER=$CXX_COMPILER
	else
		LINKER=$C_COMPILER
	fi

	if [ -n "$ANY_COMPILED_LOCAL" ]; then
		PRINT_INFO "LINKING SHARED LIB $LIB_NAME"

		ALL_FLAGS="-fPIC $C_FLAGS -shared -o $LIB_PATH ${OBJECTS[*]} ${STATIC_DEPS[*]} ${EXTRA_LIBRARY_PATH[*]} ${DYNAMIC_DEPS[*]}"

		ECHO_AND_EXEC_COMMAND "$LINKER $ALL_FLAGS"

		PRINT_INFO "COPYING SHARED LIB $LIB_NAME TO $LIBRARY_OUTPUT_DIRECTORY"
		ECHO_AND_EXEC_COMMAND "cp $LIB_PATH $LIBRARY_OUTPUT_DIRECTORY"
		touch "$TIME_FILE"
	else
		if [ ! -f "$LIBRARY_OUTPUT_DIRECTORY/$LIB_NAME".so ]; then
			PRINT_INFO "COPYING SHARED LIB $LIB_NAME TO $LIBRARY_OUTPUT_DIRECTORY"
			ECHO_AND_EXEC_COMMAND "cp $LIB_PATH $LIBRARY_OUTPUT_DIRECTORY"
		fi
	fi

	COMPILED_SHARED_LIBS[$1]=$LIB_PATH
	export COMPILED_SHARED_LIBS

}

ADD_SUBDIRECTORY() {
	local PREVIOUS_DIR
	PREVIOUS_DIR=$PWD
	cd "$1" || exit
    # shellcheck disable=SC1091
	source build.sh
	cd "$PREVIOUS_DIR" || exit
}

CHECK_CUSTOM_FILE() {

	local DEFAULT_CUSTOM_FILE="custom_functions.c"
	CUSTOM_FILE=""

	if [ -f "${DEFAULT_CUSTOM_FILE}" ]; then
		CUSTOM_FILE="${DEFAULT_CUSTOM_FILE}"
		export CUSTOM_FILE
	fi

}
