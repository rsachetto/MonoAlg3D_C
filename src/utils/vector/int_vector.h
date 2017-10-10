//
// Created by sachetto on 10/10/17.
//

#ifndef MONOALG3D_INT_VECTOR_H
#define MONOALG3D_INT_VECTOR_H

#include <stddef.h>

typedef int data_type;

/**
 * @brief A structure representing the vector object
 */
typedef struct int_vector_t {
    data_type *base; /**< Raw memory for items */
    size_t size;     /**< The number of inserted items */
    size_t capacity; /**< The number of potential items before a resize */
} int_vector;

extern int_vector *int_vector_create(size_t capacity);
extern int_vector *int_vector_clone (int_vector *v);
extern void int_vector_clear (int_vector *v);
extern void int_vector_clear_and_free_data(int_vector *v);
extern int int_vector_resize (int_vector *v, size_t capacity);
extern size_t int_vector_size (int_vector *v);
extern data_type int_vector_at (int_vector *v, size_t index);
extern int int_vector_insert (int_vector *v, data_type item, size_t index);
extern int int_vector_push_back(int_vector *v, data_type item);
extern int int_vector_remove (int_vector *v, size_t index);
extern int int_vector_remove_if (int_vector *v, int (*pred) (data_type item));
extern size_t int_vector_find (int_vector *v, data_type value);
extern size_t int_vector_find_if (int_vector *v, int (*pred) (data_type item));

#endif // MONOALG3D_INT_VECTOR_H
