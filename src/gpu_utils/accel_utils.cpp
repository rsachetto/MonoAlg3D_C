//
// Created by sachetto on 30/04/25.
//

#include "accel_utils.h"
#ifdef COMPILE_CUDA
#include "gpu_utils.h"
#endif

#ifdef COMPILE_SYCL
#include <dpct/dpct.hpp>
#include <sycl/sycl.hpp>
#endif

#include "../ode_solver/ode_solver.h"

char *get_device_info(struct ode_solver *the_ode_solver, bool linear_solver_on_gpu) {

    assert(the_ode_solver);
    char *device_info = (char *)calloc(sizeof(char), 1024);

#ifdef COMPILE_CUDA
    int device_count;
    int device = the_ode_solver->gpu_id;

    check_cuda_error(cudaGetDeviceCount(&device_count));

    if(device_count > 0) {
        if(device >= device_count) {
            sprintf(device_info, "Invalid gpu_id %d. Using gpu_id 0!\n", device);
            the_ode_solver->gpu_id = device = 0;
        }

        struct cudaDeviceProp prop;
        check_cuda_error(cudaGetDeviceProperties(&prop, the_ode_solver->gpu_id));

        if(the_ode_solver->gpu && linear_solver_on_gpu) {
            sprintf(device_info, "%d devices available, running both ODE and linear system solvers on GPU (device %d -> %s)", device_count, device, prop.name);
        } else if(the_ode_solver->gpu) {
            sprintf(device_info, "%d devices available, running only the ODE solver on GPU (device %d -> %s)", device_count, device, prop.name);
        } else {
            sprintf(device_info, "%d devices available, running only the linear system solver on GPU (device %d -> %s)", device_count, device, prop.name);
        }

        check_cuda_error(cudaSetDevice(device));
    }
#elif defined(COMPILE_SYCL)
    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();

    if(the_ode_solver->gpu && linear_solver_on_gpu) {
        sprintf(device_info, "Running both ODE and linear system solvers on GPU (%s)", q_ct1.get_device().get_info<sycl::info::device::name>().c_str());
    } else if(the_ode_solver->gpu) {
        sprintf(device_info, "Running only the ODE solver on GPU (%s)", q_ct1.get_device().get_info<sycl::info::device::name>().c_str());
    } else {
        sprintf(device_info, "Running only the linear system solver on GPU (%s)", q_ct1.get_device().get_info<sycl::info::device::name>().c_str());
    }
#endif
    return device_info;
}

void malloc_device(void **ptr, size_t n) {

#ifdef COMPILE_CUDA
    check_cuda_error(cudaMalloc(ptr, n));
#elif defined(COMPILE_SYCL)
    dpct::device_ext &dev = dpct::get_current_device();
    sycl::queue &q = dev.default_queue();

    DPCT_CHECK_ERROR(*ptr = sycl::malloc_device(n, q));
#endif
}

void free_device(void *ptr) {
#ifdef COMPILE_CUDA
    check_cuda_error(cudaFree(ptr));
#elif defined(COMPILE_SYCL)
    dpct::device_ext &dev = dpct::get_current_device();
    sycl::queue &q = dev.default_queue();
    DPCT_CHECK_ERROR(dpct::dpct_free(ptr, q));
#endif
}

void memcpy_device(void *dest, const void *src, size_t n, enum copy_direction kind) {

#ifdef COMPILE_CUDA
    enum cudaMemcpyKind dir;
    if(kind == HOST_TO_DEVICE) {
        dir = cudaMemcpyHostToDevice;
    } else if(kind == DEVICE_TO_HOST) {
        dir = cudaMemcpyDeviceToHost;
    } else if(kind == DEVICE_TO_DEVICE) {
        dir = cudaMemcpyDeviceToDevice;
    } else {
        fprintf(stderr, "Invalid copy direction\n");
        exit(EXIT_FAILURE);
    }
    check_cuda_error(cudaMemcpy(dest, src, n, dir));
#elif defined(COMPILE_SYCL)
    (void)kind;
    dpct::device_ext &dev = dpct::get_current_device();
    sycl::queue &q = dev.default_queue();
    q.memcpy(dest, src, n).wait();
#endif
}

void memcpy2d_device(void *dest, size_t pitch_dest, const void *src, size_t pitch_src, size_t w, size_t h, enum copy_direction kind) {

#ifdef COMPILE_CUDA
    enum cudaMemcpyKind dir;
    if(kind == HOST_TO_DEVICE) {
        dir = cudaMemcpyHostToDevice;
    } else if(kind == DEVICE_TO_HOST) {
        dir = cudaMemcpyDeviceToHost;
    } else if(kind == DEVICE_TO_DEVICE) {
        dir = cudaMemcpyDeviceToDevice;
    } else {
        fprintf(stderr, "Invalid copy direction\n");
        exit(EXIT_FAILURE);
    }

    check_cuda_error(cudaMemcpy2D(dest, pitch_dest, src, pitch_src, w, h, dir));
#elif defined(COMPILE_SYCL)
    (void)pitch_dest;
    (void)pitch_src;
    (void)w;
    (void)h;
    (void)kind;

    dpct::device_ext &dev = dpct::get_current_device();
    sycl::queue &q = dev.default_queue();
    q.memcpy(dest, src, w * h).wait();
#endif
}
