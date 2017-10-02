//
// Created by sachetto on 01/10/17.
//

#include "point_hash.h"
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#define INITIAL_SIZE (1024)
#define GROWTH_FACTOR (2)
#define MAX_LOAD_FACTOR (1)

/* dictionary initialization code used in both DictCreate and grow */
struct point_hash* internal_hash_create(int size)
{
    struct point_hash *d;
    int i;

    d = malloc(sizeof(*d));

    assert(d != 0);

    d->size = size;
    d->n = 0;
    d->table = malloc(sizeof(struct elt *) * d->size);

    assert(d->table != 0);

    for(i = 0; i < d->size; i++) d->table[i] = 0;

    return d;
}

struct point_hash* hash_create()
{
    return internal_hash_create(INITIAL_SIZE);
}

void hash_destroy(struct point_hash * d) {
    int i;
    struct elt *e;
    struct elt *next;

    for(i = 0; i < d->size; i++) {
        for(e = d->table[i]; e != 0; e = next) {
            next = e->next;

            //If we use a pointer we need a free, but now we dont need
            //free(e->key);

            free(e);
        }
    }

    free(d->table);
    free(d);
}

#define MULTIPLIER (97)



static unsigned long hash_function(struct point_3d k)
{

    return (unsigned long)(k.x * 18397) + (k.y * 20483) + (k.z * 29303);

}

static void grow(struct point_hash *d)  {

    struct point_hash * d2;            /* new dictionary we'll create */
    struct point_hash swap;   /* temporary structure for brain transplant */
    int i;
    struct elt *e;

    d2 = internal_hash_create(d->size * GROWTH_FACTOR);

    for(i = 0; i < d->size; i++) {
        for(e = d->table[i]; e != 0; e = e->next) {
            /* note: this recopies everything */
            /* a more efficient implementation would
             * patch out the strdups inside DictInsert
             * to avoid this problem */
            hash_insert(d2, e->key, e->value);
        }
    }

    /* the hideous part */
    /* We'll swap the guts of d and d2 */
    /* then call DictDestroy on d2 */
    swap = *d;
    *d = *d2;
    *d2 = swap;

    hash_destroy(d2);
}

/* insert a new key-value pair into an existing dictionary */
void hash_insert(struct point_hash* d, struct point_3d key, int value)
{
    struct elt *e;
    unsigned long h;

    e = malloc(sizeof(*e));

    assert(e);

    e->key.x = key.x;
    e->key.y = key.y;
    e->key.z = key.z;

    e->value = value;

    h = hash_function(key) % d->size;

    e->next = d->table[h];
    d->table[h] = e;

    d->n++;

    /* grow table if there is not enough room */
    if(d->n >= d->size * MAX_LOAD_FACTOR) {
        grow(d);
    }
}

bool point_equals(struct point_3d a, struct point_3d b) {
    return (a.x == b.x) && (a.y == b.y) && (a.z == b.z);
}

/* return the most recently inserted value associated with a key */
/* or 0 if no matching key is present */
int hash_search(struct point_hash* d, struct point_3d key) {
    struct elt *e;

    for(e = d->table[hash_function(key) % d->size]; e != 0; e = e->next) {
        if(point_equals(e->key, key)) {
            /* got it */
            return e->value;
        }
    }

    return 0;
}

/* delete the most recently inserted record with the given key */
/* if there is no such record, has no effect */
void hash_delete(struct point_hash *d, struct point_3d key)
{
    struct elt **prev;          /* what to change when elt is deleted */
    struct elt *e;              /* what to delete */

    for(prev = &(d->table[hash_function(key) % d->size]);
        *prev != 0;
        prev = &((*prev)->next)) {
        if(point_equals((*prev)->key, key)) {
            /* got it */
            e = *prev;
            *prev = e->next;

            //free(e->value);
            free(e);

            return;
        }
    }
}