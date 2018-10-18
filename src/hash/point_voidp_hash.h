//
// Created by sachetto on 01/10/17.
//

#ifndef MONOALG3D_POINT_VOID_HASH_H
#define MONOALG3D_POINT_VOID_HASH_H

#include "hash_common.h"

struct elt {
    struct elt *next;
    struct point_3d key;
    void *value;
};

struct point_voidp_hash {
    int size;           /* size of the pointer table */
    int n;              /* number of elements stored */
    struct elt **table;
};

/* create a new empty dictionary */
struct point_voidp_hash* point_voidp_hash_create();

/* destroy a dictionary */
void point_voidp_hash_destroy(struct point_voidp_hash*);

/* insert a new key-value pair into an existing dictionary */
void point_voidp_hash_insert(struct point_voidp_hash*, struct point_3d key, void* value);

/* return the most recently inserted value associated with a key */
/* or 0 if no matching key is present */
void* point_voidp_hash_search(struct point_voidp_hash*, struct point_3d key);

/* delete the most recently inserted record with the given key */
/* if there is no such record, has no effect */
void point_voidp_hash_delete(struct point_voidp_hash*, struct point_3d key);

#endif //MONOALG3D_HASH_H
