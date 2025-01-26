#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include "gpu_utils.h"

void kernel_gpu_vec_div_vec(real *vec1, real *vec2, real *vec3, size_t n, const sycl::nd_item<3> &item_ct1) {

    unsigned int i = item_ct1.get_group(2) * item_ct1.get_local_range(2) + item_ct1.get_local_id(2);

    if(i < n) {
        vec3[i] = vec1[i]/vec2[i];
    }
}

extern "C" void gpu_vec_div_vec(real *vec1, real *vec2, real *res, size_t n) {
    const int GRID  = (n + BLOCK_SIZE - 1)/BLOCK_SIZE;
    dpct::get_in_order_queue().parallel_for(sycl::nd_range<3>(sycl::range<3>(1, 1, GRID) * sycl::range<3>(1, 1, BLOCK_SIZE), sycl::range<3>(1, 1, BLOCK_SIZE)),
                                            [=](sycl::nd_item<3> item_ct1) {
                                                kernel_gpu_vec_div_vec(vec1, vec2, res, n, item_ct1);
                                            });
    dpct::get_current_device().queues_wait_and_throw();
}