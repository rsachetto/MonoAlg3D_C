//
// Created by sachetto on 10/10/17.
//



#include <stdlib.h>
#include <string.h>
#include "int_vector.h"

static int auto_capacity(int_vector *v, size_t *new_capacity);

/**
 * @brief Creates a new vector object
 * @param[in] item_size The size of each element in bytes
 * @param[in] capacity Default number of items allocated to the vector
 * @return A pointer to the new vector object, or NULL on failure
 * @remarks If capacity is 0, an empty vector will be created
 */
extern int_vector *int_vector_create(size_t capacity)
{
    int_vector *v = malloc(sizeof *v);

    if (v != NULL) {
        v->base = NULL;
        v->size = 0;
        v->capacity = capacity;

        if (capacity > 0) {
            /* Allocate the default capacity */
            v->base = malloc(capacity * sizeof(data_type));

            if (v->base == NULL) {
                /* Clean up rather than leaving a zombie */
                free(v);
                v = NULL;
            }
        }
    }

    return v;
}

/**
 * @brief Creates an exact copy of vector object
 * @param[in] v The vector being copied
 * @return A pointer to the new vector object, or NULL on failure
 */
extern int_vector *int_vector_clone(int_vector *v)
{
    int_vector *p = int_vector_create(v->capacity);

    if (p != NULL) {
        /* Copy the parts that int_vector_create() didn't */
        memcpy(p->base, v->base, v->size * sizeof(data_type));
        p->size = v->size;
    }

    return p;
}

/**
 * @brief Clears a vector object created by vector_new() to an empty state
 * @param[in] v A pointer to the vector being cleared
 */
extern void int_vector_clear(int_vector *v)
{
    v->size = 0;
}

/**
 * @brief Clears a vector object created by vector_new() to an empty state
 * @param[in] v A pointer to the vector being destroyed
 */
extern void int_vector_clear_and_free_data(int_vector *v)
{
    v->size = 0;
    v->capacity = 0;
    free(v->base);
}

/**
 * @brief Resizes a vector object to the specified capacity
 * @param[in] v A pointer to the vector being resized
 * @param[in] capacity The new capacity
 * @return True/non-zero if the resize succeeded, false/zero otherwise
 * @remarks 1) The new capacity may be larger or smaller
 *          2) No-op if the capacity is unchanged, with a success return
 *          3) Resizing may invalidate pointers into the vector
 */
extern int int_vector_resize(int_vector *v, size_t capacity)
{
    if (capacity == 0)
        return 0;

    if (capacity != v->capacity) {
        void *temp;

        v->capacity = capacity;

        /*
            If the int_vector is empty, realloc() depends on v->base
            being initialized to NULL
        */
        temp = realloc(v->base, v->capacity * sizeof(data_type));

        if (temp == NULL)
            return 0;

        v->base = temp;
    }

    return 1;
}

/**
 * @brief Returns the size of the specified vector
 * @param[in] v A pointer to the vector
 * @return The number of items inserted into the vector
 */
extern size_t int_vector_size(int_vector *v)
{
    return v->size;
}

/**
 * @brief Returns the item stored at the specified index of a vector
 * @param[in] v A pointer to the vector
 * @param[in] index The index of the requested item
 * @return A generic pointer to the item
 */
extern data_type int_vector_at(int_vector *v, size_t index)
{
    return v->base[index];
}

/**
 * @brief Inserts a single item in a vector at the specified index
 * @param[in] v A pointer to the vector being inserted into
 * @param[in] item A pointer to the item being appended
 * @param[in] index The index where the item will be inserted
 * @return True/non-zero if the insert succeeded, false/zero otherwise
 * @remarks 1) The vector may be resized to make room
 *          2) The item pointed to is copied by value into the vector
 *          3) The size of the vector will increase by 1 on success
 *          4) The capacity of the vector may increase by more than 1
 *          5) All items from index to v->size may be shifted to make room
 */
extern int int_vector_insert(int_vector *v, data_type item, size_t index)
{
    int *src, *dst;

    if (index > v->size)
        return 0;

    if (v->size == v->capacity) {
        /* Resize to the next auto-growth amount */
        size_t new_capacity;

        if (!auto_capacity(v, &new_capacity) ||
            !int_vector_resize(v, new_capacity))
        {
            return 0;
        }

        v->capacity = new_capacity;
    }

    src = &v->base[index + 1];
    dst = &v->base[index];

    /* Make room for a new item */
    memmove(src, dst, (v->size - index) * sizeof(data_type));

    /* Insert the new item */
    memcpy(dst, &item, sizeof(data_type));
    ++v->size;

    return 1;
}

/**
 * @brief Inserts a single item in the end of the vector
 * @param[in] v A pointer to the vector being inserted into
 * @param[in] item A pointer to the item being appended
 * @return True/non-zero if the insert succeeded, false/zero otherwise
 *
 */
extern int int_vector_push_back(int_vector *v, data_type item)
{

    if (v->size == v->capacity) {
        /* Resize to the next auto-growth amount */
        size_t new_capacity;

        if (!auto_capacity(v, &new_capacity) ||
            !int_vector_resize(v, new_capacity))
        {
            return 0;
        }

        v->capacity = new_capacity;
    }

    v->base[v->size] = item;

    ++v->size;

    return 1;
}



/**
 * @brief Removes the specified item within the a vector
 * @param[in] v A pointer to the vector
 * @param[in] index The index of the item
 * @return True/non-zero if the value was found and removed, false/zero otherwise
 * @remarks All items following the found value may be shifted to fill in the hole
 */
extern int int_vector_remove(int_vector *v, size_t index)
{
    if (index >= v->size)
        return 0;
    else if (index == v->size - 1) {
        /* Special case: no copy when removing the last item */
        --v->size;
    }
    else {
        int *dst =  &v->base[index];
        int *src =  &v->base[index + 1];

        /* Fill in the vacated slot */
        memmove(dst, src, (v->size - index - 1) * sizeof(data_type));
        --v->size;
    }

    return 1;
}

/**
 * @brief Removes the specified item within the a vector
 * @param[in] v A pointer to the vector
 * @param[in] pred A predicate for determining the first matching item
 * @return True/non-zero if the value was found and removed, false/zero otherwise
 * @remarks All items following the found value may be shifted to fill in the hole
 */
extern int int_vector_remove_if(int_vector *v, int (*pred)(data_type item))
{
    size_t index = int_vector_find_if(v, pred);

    if (index != -1)
        return int_vector_remove(v, index);

    return 0;
}

/**
 * @brief Searches for the specified value within the a vector
 * @param[in] v A pointer to the vector
 * @param[in] item the value being searched for
 * @return The index of the found value, or (size_t)-1 if not found
 *  */
extern size_t int_vector_find(int_vector *v, data_type  value)
{
    size_t i;

    for (i = 0; i < v->size; i++) {
        if ( v->base[i] == value)
            return i;
    }

    return -1;
}

/**
 * @brief Searches for the specified value within the a vector
 * @param[in] v A pointer to the vector
 * @param[in] pred A predicate for determining the first matching item
 * @return The index of the found value, or (size_t)-1 if not found
 */
extern size_t int_vector_find_if(int_vector *v, int (*pred)(data_type item))
{
    size_t i;

    for (i = 0; i < v->size; i++) {
        if (pred(int_vector_at(v, i)))
            return i;
    }

    return -1;
}

/**
 * @brief Calculates the auto-growth of a vector
 * @param[in] v A pointer to the vector
 * @param[out] new_capacity The calculated new capacity
 * @return True/non-zero if the vector can grow, false/zero otherwise
 */
static int auto_capacity(int_vector *v, size_t *new_capacity)
{
    if (v->capacity == -1)
        return 0;
    else if (v->capacity == 0) {
        /* Default to something reasonable for an empty int_vector */
        *new_capacity = 4;
    }
    else {
        /* Try to grow by 50% of the current capacity */
        *new_capacity = v->capacity * 1.5;

        if (*new_capacity < v->capacity) {
            /* Max out the capacity if 50% growth overflows (wraps) */
            *new_capacity = -1;
        }
    }

    return 1;
}