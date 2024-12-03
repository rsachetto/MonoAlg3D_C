//
// Created by sachetto on 28/11/24.
//

#ifndef ACCEL_UTILS_H
#define ACCEL_UTILS_H

#include <stddef.h>

typedef enum {
  HOST_TO_DEVICE,
  DEVICE_TO_HOST,
} copy_direction;

#ifdef __cplusplus
extern "C" {
#endif
  void malloc_device(void **ptr, size_t n);
  void free_device(void *ptr);
  void memcpy_device(void *dest, const void *src, size_t n, copy_direction kind);
  void create_sparse_handle(void *handle);
  void create_blas_handle(void *handle);
#ifdef __cplusplus
}
#endif
#endif //ACCEL_UTILS_H
