//
// Created by sachetto on 30/04/25.
//

#ifndef ACCEL_UTILS_H
#define ACCEL_UTILS_H

#include <stddef.h>
#include <stdint.h>

enum copy_direction { HOST_TO_DEVICE, DEVICE_TO_HOST, DEVICE_TO_DEVICE };

#ifdef __cplusplus
extern "C" {
#endif
void malloc_device(void **ptr, size_t n);
void free_device(void *ptr);
void memcpy_device(void *dest, const void *src, size_t n, enum copy_direction kind);
void memcpy2d_device(void *dest, size_t pitch_dest, const void *src, size_t pitch_src, size_t w, size_t h, enum copy_direction kind);
#ifdef __cplusplus
}
#endif

#endif // ACCEL_UTILS_H
