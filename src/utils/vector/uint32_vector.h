//
// Created by sachetto on 10/10/17.
//

#ifndef MONOALG3D_U32_vector_H
#define MONOALG3D_U32_vector_H

#include <stddef.h>

typedef u_int32_t data_type_32;

/**
 * @brief A structure representing the vector object
 */
typedef struct uint32_vector_t {
    data_type_32 *base; /**< Raw memory for items */
    size_t size;     /**< The number of inserted items */
    size_t capacity; /**< The number of potential items before a resize */
} uint32_vector;

 uint32_vector *uint32_vector_create(size_t capacity);
 uint32_vector *uint32_vector_clone (uint32_vector *v);
 void uint32_vector_clear (uint32_vector *v);
 void uint32_vector_clear_and_free_data(uint32_vector *v);
 int uint32_vector_resize (uint32_vector *v, size_t capacity);
 size_t uint32_vector_size (uint32_vector *v);
 data_type_32 uint32_vector_at (uint32_vector *v, size_t index);
 int uint32_vector_insert (uint32_vector *v, data_type_32 item, size_t index);
 int uint32_vector_push_back(uint32_vector *v, data_type_32 item);
 data_type_32 uint32_vector_pop_back(uint32_vector *v);
 int uint32_vector_remove (uint32_vector *v, size_t index);
 int uint32_vector_remove_if (uint32_vector *v, int (*pred) (data_type_32 item));
 size_t uint32_vector_find (uint32_vector *v, data_type_32 value);
 size_t uint32_vector_find_if (uint32_vector *v, int (*pred) (data_type_32 item));

#endif // MONOALG3D_U32_vector_H
