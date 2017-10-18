//
// Created by sachetto on 01/10/17.
//

#ifndef MONOALG3D_STRING_HASH_H
#define MONOALG3D_STRING_HASH_H

#define STRING_HASH_PRINT_KEY_VALUE(d)                                                                          \
    for (int i = 0; i < (d)->size; i++) {                                                                                   \
        for (struct string_elt *e = (d)->table[i % (d)->size]; e != 0; e = e->next) {                                      \
            printf("%s = %s\n", e->key, e->value);                                                                                             \
        }                                                                                                              \
    }

#define STRING_HASH_PRINT_KEY_VALUE_LOG(d)                                                                          \
    for (int i = 0; i < (d)->size; i++) {                                                                                   \
        for (struct string_elt *e = (d)->table[i % (d)->size]; e != 0; e = e->next) {                                      \
            print_to_stdout_and_file("%s = %s\n", e->key, e->value);                                                                                             \
        }                                                                                                              \
    }

struct string_elt {
    struct string_elt *next;
    char *key;
    char *value;
};


struct string_hash {
    int size; /* size of the pointer table */
    int n;    /* number of elements stored */
    struct string_elt **table;
};

/* create a new empty dictionary */
struct string_hash *string_hash_create ();

/* destroy a dictionary */
void string_hash_destroy (struct string_hash *);

/* insert a new key-value pair into an existing dictionary */
void string_hash_insert (struct string_hash *, const char *key, const char *value);

/* return the most recently inserted value associated with a key */
/* or NULL if no matching key is present */
char *string_hash_search (struct string_hash *, const char *key);

/* delete the most recently inserted record with the given key */
/* if there is no such record, has no effect */
void string_hash_delete (struct string_hash *, const char *key);

void string_hash_insert_or_overwrite(struct string_hash *d, const char *key, const char *value);

#endif // MONOALG3D_STRING_HASH_H
