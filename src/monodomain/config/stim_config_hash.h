//
// Created by sachetto on 01/10/17.
//

#ifndef MONOALG3D_STIM_CONFIG_HASH_H
#define MONOALG3D_STIM_CONFIG_HASH_H

#include "stim_config.h"

#define STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE(d, fn)                                                                          \
    for (int i = 0; i < (d)->size; i++) {                                                                                   \
        for (struct stim_config_elt *e = (d)->table[i % (d)->size]; e != 0; e = e->next) {                                      \
            fn (e->value);                                                                                             \
        }                                                                                                              \
    }

#define STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE_AND_KEY(d, fn)                                                                          \
    for (int i = 0; i < (d)->size; i++) {                                                                                   \
        for (struct stim_config_elt *e = (d)->table[i % (d)->size]; e != 0; e = e->next) {                                      \
            fn (e->value, e->key);                                                                                             \
        }                                                                                                              \
    }

struct stim_config_elt {
    struct stim_config_elt *next;
    char *key;
    struct stim_config *value;
};

struct stim_config_hash {
    size_t size; /* size of the pointer table */
    size_t n;    /* number of elements stored */
    struct stim_config_elt **table;
};


/* create a new empty dictionary */
struct stim_config_hash *stim_config_hash_create();

/* destroy a dictionary */
void stim_config_hash_destroy (struct stim_config_hash *);

/* insert a new key-value pair into an existing dictionary */
void stim_config_hash_insert (struct stim_config_hash *, const char *key, struct stim_config *value);

/* return the most recently inserted value associated with a key */
/* or NULL if no matching key is present */
struct stim_config *stim_config_hash_search (const struct stim_config_hash *, const char *key);

/* delete the most recently inserted record with the given key */
/* if there is no such record, has no effect */
void stim_config_hash_delete (struct stim_config_hash *, const char *key);

#endif // MONOALG3D_STIM_CONFIG_HASH_H
