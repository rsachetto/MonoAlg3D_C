SET(LIB_TYPE SHARED) #set the lib type

set(MODEL_FILE_CPU "bondarenko_2004_CPU.c")
set(MODEL_FILE_GPU "bondarenko_2004_GPU.cu")

set(MODEL_FILE_CPU "ten_tusscher_2006_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2006_GPU.cu")

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(bondarenko_2004 ${LIB_TYPE} model_common.h ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
    CUDA_ADD_LIBRARY(ten_tusscher_2006 ${LIB_TYPE} model_common.h ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
    ADD_LIBRARY(bondarenko_2004 ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
    ADD_LIBRARY(ten_tusscher_2006 ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
endif(CUDA_FOUND)




