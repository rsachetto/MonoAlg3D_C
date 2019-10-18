#!/bin/bash

FIND_CUDA () {

	#CUDA RELATED VARIABLES
	CUDA_LIBRARY_PATH=""
	CUDA_INCLUDE_PATH=""
	NVCC=""
	CUDA_FOUND=""

	if [ -z "$CUDA_LIBRARY_PATH" ]; then
		CUDA_LIBRARY_PATH=$(dirname "$(ldconfig -p | grep libcudart | awk '{print $4}' | head -n 1 | head -c -5)" 2> /dev/null)
	fi

	if [ -z "$CUDA_INCLUDE_PATH" ]; then
		CUDA_INCLUDE_PATH=${CUDA_LIBRARY_PATH/lib64/include}
	fi

	if [ -z "$NVCC" ]; then
		NVCC=$(command -v nvcc | xargs) #in some systems command stops the script
	fi

	if [ ! -d "$CUDA_LIBRARY_PATH" ] || [ ! -d "$CUDA_INCLUDE_PATH" ] || [ ! -f "$NVCC" ]; then

	if [ ! -d "$CUDA_LIBRARY_PATH" ] ||  [ ! -f "$CUDA_LIBRARY_PATH/libcudart.so" ]; then
	
		if [ -z "$CUDA_LIBRARY_PATH" ]; then
			PRINT_WARN "CUDA_LIBRARY_PATH variable is empty. CUDA code will not be compiled!"
		else
			PRINT_WARN "CUDA_LIBRARY_PATH set as $CUDA_LIBRARY_PATH is not a valid CUDA library path. CUDA code will not be compiled!"
		fi
		
		PRINT_WARN "You can change the CUDA_LIBRARY_PATH variable manually in the file $OPTIONS_FILE"
		echo
	fi

	if [ ! -d "$CUDA_INCLUDE_PATH" ] ||  [ ! -f "$CUDA_INCLUDE_PATH/cuda_runtime.h" ]; then
	
		if [ -z "$CUDA_INCLUDE_PATH" ]; then
			PRINT_WARN "CUDA_INCLUDE_PATH variable is empty. CUDA code will not be compiled!"
		else
			PRINT_WARN "CUDA_INCLUDE_PATH set as $CUDA_INCLUDE_PATH is not a valid CUDA include path. CUDA code will not be compiled!"
		fi
		
		PRINT_WARN "You can change the CUDA_INCLUDE_PATH variable manually in the file $OPTIONS_FILE"
		echo
	fi

	if [ ! -f "$NVCC" ]; then
		PRINT_WARN "CUDA compiler set as $NVCC not found. CUDA code will not be compiled!"
		PRINT_WARN "You can change the NVCC variable manually in the file $OPTIONS_FILE"
	fi
	else
		CUDA_FOUND='y'
		CUDA_LIBRARIES="cudart"
	fi
}

FIND_CRITERION() {

	PKG_CONFIG=$(command -v pkg-config)	

	if [ -z "$PKG_CONFIG" ]; then
		CRITERION_FOUND=''
	else
		local EXISTS
		$PKG_CONFIG --exists criterion
		EXISTS=$?
		if [ "$EXISTS" == 0 ]; then
			CRITERION_FOUND='y'
			CRITERIONS_LIBRARY_PATH=$($PKG_CONFIG --libs-only-L criterion)
			CRITERION_INCLUDE_PATH=$($PKG_CONFIG --cflags-only-I criterion)
			LIBS=$($PKG_CONFIG --libs-only-l criterion)
			
			if [ -n "$LIBS" ]; then
				for i in $LIBS; do
					i=${i:2}
					CRITERION_LIBRARIES="$CRITERION_LIBRARIES $i"
				done			
			fi

		fi	
	fi

}

FIND_MPI () {
	#TODO: command to find MPI
	MPI_FOUND='y'
	MPI_LIBRARIES=mpi
	MPI_LIBRARY_PATH=/usr/lib/openmpi/
	MPI_INCLUDE_PATH=''
}