//
// Created by sachetto on 01/10/17.
//

#include "stim_config_hash.h"
#include <assert.h>
#include <string.h>

#include "../../hash/hash_common.h"


/* dictionary initialization code used in both create_hash and grow */
struct stim_config_hash *internal_stim_config_hash_create(size_t size) {
    struct stim_config_hash *d;
    int i;

    d = malloc (sizeof (*d));

    assert (d != 0);

    d->size = size;
    d->n = 0;
    d->table = malloc (sizeof (struct stim_config_elt *) * d->size);

    assert (d->table != 0);

    for (i = 0; i < d->size; i++)
        d->table[i] = NULL;

    return d;
}

struct stim_config_hash *stim_config_hash_create () {
    return internal_stim_config_hash_create (INITIAL_SIZE);
}

void stim_config_hash_destroy (struct stim_config_hash *d) {
    int i;
    struct stim_config_elt *e;
    struct stim_config_elt *next;

    for (i = 0; i < d->size; i++) {
        for (e = d->table[i]; e != 0; e = next) {
            next = e->next;
            free(e->key);
            //we dont free the value as it does not belong to us!
            free (e);
        }
    }

    free (d->table);
    free (d);
}

static unsigned long stim_config_hash_function (const char *k) {

    unsigned long hash = 5381;
    int c;

    while ((c = *k++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

static void grow (struct stim_config_hash *d) {

    struct stim_config_hash *d2;  /* new dictionary we'll create */
    struct stim_config_hash swap; /* temporary structure for brain transplant */
    int i;
    struct stim_config_elt *e;

    d2 = internal_stim_config_hash_create (d->size * GROWTH_FACTOR);

    for (i = 0; i < d->size; i++) {
        for (e = d->table[i]; e != 0; e = e->next) {
            /* note: this recopies everything */
            stim_config_hash_insert (d2, e->key, e->value);
        }
    }

    /* the hideous part */
    /* We'll swap the guts of d and d2 */
    swap = *d;
    *d = *d2;
    *d2 = swap;

    stim_config_hash_destroy (d2);
}

/* insert a new key-value pair into an existing dictionary */
void stim_config_hash_insert (struct stim_config_hash *d, const char *key,  struct stim_config *value) {
    struct stim_config_elt *e;
    unsigned long h;

    e = malloc (sizeof (*e));

    assert (e);

    e->key = strdup (key);

    e->value = value;

    h = stim_config_hash_function (key) % d->size;

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
struct stim_config *stim_config_hash_search (const struct stim_config_hash *d, const char *key) {
    struct stim_config_elt *e;

    for (e = d->table[stim_config_hash_function (key) % d->size]; e != 0; e = e->next) {
        if (strcmp (e->key, key) == 0) {
            /* got it */
            return e->value;
        }
    }

    return NULL;
}

/* delete the most recently inserted record with the given key */
/* if there is no such record, has no effect */
void stim_config_hash_delete (struct stim_config_hash *d, const char *key) {
    struct stim_config_elt **prev; /* what to change when elt is deleted */
    struct stim_config_elt *e;     /* what to delete */

    for (prev = &(d->table[stim_config_hash_function (key) % d->size]); *prev != 0; prev = &((*prev)->next)) {
        if (strcmp ((*prev)->key, key) == 0) {
            /* got it */
            e = *prev;
            *prev = e->next;

            free (e->key);
            free (e);

            return;
        }
    }
}