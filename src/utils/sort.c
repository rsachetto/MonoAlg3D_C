//
// Created by sachetto on 01/10/17.
//

#include "utils.h"
#define SWAP(x, y, T) do { T SWAP = (x); (x) = (y); (y) = SWAP; } while (0)

int partition(double **a, int column, int first, int last) {
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
    double pivot = a[first][column];

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
            SWAP(a[i], a[loc], double*);
        }
    }

    // i greater than loc, switch first w/ loc to sort
    SWAP(a[first], a[loc], double*);

    return loc; // return current loc
}

void quicksort(double **a, int column, int first, int last) {
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


void sort_vector(double **a, int length) {

    quicksort(a, 0, 0, length - 1); // perform quicksort algorithm ordering using the X coordinate

    int j = 0;
    while ( j < length ) {
        double xnum = a[j][0];
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
        double xnum = a[k][0];
        double ynum = a[k][1];
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

