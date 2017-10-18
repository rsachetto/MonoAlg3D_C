//
// Created by sachetto on 01/10/17.
//

#include "string_hash.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "hash_common.h"

/* dictionary initialization code used in both create_hash and grow */
struct string_hash *internal_string_hash_create (int size) {
    struct string_hash *d;
    int i;

    d = malloc (sizeof (*d));

    assert (d != 0);

    d->size = size;
    d->n = 0;
    d->table = malloc (sizeof (struct string_elt *) * d->size);

    assert (d->table != 0);

    for (i = 0; i < d->size; i++)
        d->table[i] = NULL;

    return d;
}

struct string_hash *string_hash_create () {
    return internal_string_hash_create (INITIAL_SIZE);
}

void string_hash_destroy (struct string_hash *d) {
    int i;
    struct string_elt *e;
    struct string_elt *next;

    for (i = 0; i < d->size; i++) {
        for (e = d->table[i]; e != 0; e = next) {
            next = e->next;

            free(e->key);
            free(e->value);
            free (e);
        }
    }

    free (d->table);
    free (d);
}

static unsigned long string_hash_function (const char *k) {

    unsigned long hash = 5381;
    int c;

    while ((c = *k++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

static void grow (struct string_hash *d) {

    struct string_hash *d2;  /* new dictionary we'll create */
    struct string_hash swap; /* temporary structure for brain transplant */
    int i;
    struct string_elt *e;

    d2 = internal_string_hash_create (d->size * GROWTH_FACTOR);

    for (i = 0; i < d->size; i++) {
        for (e = d->table[i]; e != 0; e = e->next) {
            /* note: this recopies everything */
            string_hash_insert (d2, e->key, e->value);
        }
    }

    /* the hideous part */
    /* We'll swap the guts of d and d2 */
    swap = *d;
    *d = *d2;
    *d2 = swap;

    string_hash_destroy (d2);
}

void string_hash_insert_or_overwrite(struct string_hash *d, const char *key, const char *value) {

    char* new_key = string_hash_search(d, key);

    //Key found, delete it first
    if(new_key != NULL) {
        string_hash_delete(d, key);
    }

    string_hash_insert(d, key, value);


}


/* insert a new key-value pair into an existing dictionary */
void string_hash_insert (struct string_hash *d, const char *key, const char *value) {
    struct string_elt *e;
    unsigned long h;

    e = malloc (sizeof (*e));

    assert (e);

    e->key = strdup (key);

    e->value = strdup (value);

    h = string_hash_function (key) % d->size;

    e->next = d->table[h];
    d->table[h] = e;

    d->n++;

    /* grow table if there is not enough room */
    if (d->n >= d->size * MAX_LOAD_FACTOR) {
        grow (d);
    }
}

/* return the most recently inserted value associated with a key */
/* or NULL if no matching key is present */
char *string_hash_search (struct string_hash *d, const char *key) {
    struct string_elt *e;

    for (e = d->table[string_hash_function (key) % d->size]; e != 0; e = e->next) {
        if (strcmp (e->key, key) == 0) {
            /* got it */
            return strdup (e->value);
        }
    }

    return NULL;
}

/* delete the most recently inserted record with the given key */
/* if there is no such record, has no effect */
void string_hash_delete (struct string_hash *d, const char *key) {
    struct string_elt **prev; /* what to change when elt is deleted */
    struct string_elt *e;     /* what to delete */

    for (prev = &(d->table[string_hash_function (key) % d->size]); *prev != 0; prev = &((*prev)->next)) {
        if (strcmp ((*prev)->key, key) == 0) {
            /* got it */
            e = *prev;
            *prev = e->next;

            free (e->key);
            free (e->value);
            free (e);

            return;
        }
    }
}