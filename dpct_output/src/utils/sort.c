//
// Created by sachetto on 01/10/17.
//

#include "utils.h"
#define SWAP(x, y, T) do { T SWAP = (x); (x) = (y); (y) = SWAP; } while (0)

int partition(real_cpu **a, int column, int first, int last) {
    //
    // partitions array from first to last index w/ a[first] as pivot.
    // pre: a is an array where elements a[first] through a[last] have values.
    // first and last are indices within the range of a
    // post: a[first..last] has been partition where a[first..(loc-1)]
    // contains elements <= pivot; a[loc] = pivot; and a[(loc+1)..last]
    // contains elements >= pivot.
    // loc is location of pivot in partitioned array
    //

    int  i = first;
    int loc = last + 1;
    real_cpu pivot = a[first][column];

    // increase i until greater than loc
    while (i < loc) {
        // increase i while decreasing loc
        do {
            ++i;
        }
        while ((a[i][column] < pivot) && (i < last));

        do {
            --loc;
        }
        while (a[loc][column] > pivot);

        // swap a[i] and a[loc] for correct order relative to pivot
        if (i < loc) {
            SWAP(a[i], a[loc], real_cpu*);
        }
    }

    // i greater than loc, switch first w/ loc to sort
    SWAP(a[first], a[loc], real_cpu*);

    return loc; // return current loc
}

void quicksort(real_cpu **a, int column, int first, int last) {
    //
    // quicksort algorithm to sort an array of elements
    // pre: a is an array where elements a[first] through a[last] have values.
    // first and last are indices within the range of a
    // post: a[first..last] has been sorted
    //

    int loc;

    // partition everything, partition left, partition right
    // continue until only one element, now in correct index
    if (first < last) {
        loc = partition(a, column, first, last);
        quicksort(a, column, first, loc - 1);
        quicksort(a, column, loc + 1, last);
    }
}


void sort_vector(real_cpu **a, int length) {

    quicksort(a, 0, 0, length - 1); // perform quicksort algorithm ordering using the X coordinate

    int j = 0;
    while ( j < length ) {
        real_cpu xnum = a[j][0];
        int firstx = j;
        while (a[j][0] == xnum) {
            j++;
            if (j == length)
                break;
        }
        int lastx = j - 1;
        quicksort(a, 1, firstx, lastx);
    }
    int k = 0;
    while ( k < length ) {
        real_cpu xnum = a[k][0];
        real_cpu ynum = a[k][1];
        int firstxy = k;
        while (a[k][0] == xnum && a[k][1] == ynum) {
            k++;
            if (k == length) {
                break;
            }
        }
        int lastxy = k - 1;
        quicksort(a, 2, firstxy, lastxy);
    }
}

int partition_by_distance (real_cpu *dist_array, uint32_t *indexes, int first, int last) {

    int  i = first;
    int loc = last + 1;
    real_cpu pivot = dist_array[first];

    // increase i until greater than loc
    while (i < loc) {
        // increase i while decreasing loc
        do {
            ++i;
        }
        while ((dist_array[i] < pivot) && (i < last));

        do {
            --loc;
        }
        while (dist_array[loc] > pivot);

        // swap a[i] and a[loc] for correct order relative to pivot
        if (i < loc) {
            SWAP(dist_array[i], dist_array[loc], real_cpu);
            SWAP(indexes[i], indexes[loc], uint32_t);
        }
    }

    // i greater than loc, switch first w/ loc to sort
    SWAP(dist_array[first], dist_array[loc], real_cpu);
    SWAP(indexes[first], indexes[loc], uint32_t);

    return loc; // return current loc
}

void quicksort_by_distance(real_cpu *dist_array, uint32_t *indexes, int first, int last) {

    int loc;

    if (first < last) {
        loc = partition_by_distance(dist_array, indexes, first, last);
        quicksort_by_distance(dist_array, indexes, first, loc - 1);
        quicksort_by_distance(dist_array, indexes, loc + 1, last);
    }
}

void sort_vector_by_distance (real_cpu *dist_array, uint32_t *indexes, int length) {

    quicksort_by_distance(dist_array, indexes, 0, length - 1);
}

