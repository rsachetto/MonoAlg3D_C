#include "gpu_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static const char * cudaGetErrorEnum(int code, const char *api) {
	if(strcmp(api, "cuda") == 0) {
		cudaError_t error = (cudaError_t) code;
		return  cudaGetErrorName(error);
	}
	else if(strcmp(api, "cublas") == 0) {
#ifdef CUBLAS_API_H_
		cublasStatus_t error = (cublasStatus_t) code;

		switch (error) {
			case CUBLAS_STATUS_SUCCESS:
				return "CUBLAS_STATUS_SUCCESS";

			case CUBLAS_STATUS_NOT_INITIALIZED:
				return "CUBLAS_STATUS_NOT_INITIALIZED";

			case CUBLAS_STATUS_ALLOC_FAILED:
				return "CUBLAS_STATUS_ALLOC_FAILED";

			case CUBLAS_STATUS_INVALID_VALUE:
				return "CUBLAS_STATUS_INVALID_VALUE";

			case CUBLAS_STATUS_ARCH_MISMATCH:
				return "CUBLAS_STATUS_ARCH_MISMATCH";

			case CUBLAS_STATUS_MAPPING_ERROR:
				return "CUBLAS_STATUS_MAPPING_ERROR";

			case CUBLAS_STATUS_EXECUTION_FAILED:
				return "CUBLAS_STATUS_EXECUTION_FAILED";

			case CUBLAS_STATUS_INTERNAL_ERROR:
				return "CUBLAS_STATUS_INTERNAL_ERROR";

			case CUBLAS_STATUS_NOT_SUPPORTED:
				return "CUBLAS_STATUS_NOT_SUPPORTED";

			case CUBLAS_STATUS_LICENSE_ERROR:
				return "CUBLAS_STATUS_LICENSE_ERROR";
		}

		return "<unknown>";
#else
		return "<unknown>";
#endif

	}

	return "<unknown>";
}

void cuda_assert(int code, char const *const func, const char *const file, int const line, const char *api) {

	#ifdef DEBUG_INFO
	if (code) {
		fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line, code, cudaGetErrorEnum(code, api), func);
		exit(code);
	}
	#endif
}




