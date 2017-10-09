SET(LIB_TYPE SHARED) #set the lib type

set(MODEL_FILE_CPU "bondarenko_2004_CPU.c")
set(MODEL_FILE_GPU "bondarenko_2004_GPU.cu")

ADD_LIBRARY(bondarenko_2004_cpu ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
CUDA_ADD_LIBRARY(bondarenko_2004_gpu ${LIB_TYPE} ${model_common.h} ${MODEL_FILE_GPU})

set(MODEL_FILE_CPU "ten_tusscher_2006_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2006_GPU.cu")

ADD_LIBRARY(ten_tusscher_2006_cpu ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
CUDA_ADD_LIBRARY(ten_tusscher_2006_gpu ${LIB_TYPE} ${model_common.h} ${MODEL_FILE_GPU})