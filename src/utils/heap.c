#include <stdio.h>
#include <stdlib.h>
#include "heap.h"

static void heapify(struct point_distance_heap* h, int index) {
    int left = index * 2 + 1;
    int right = index * 2 + 2;
    int min = index;
 
    if (left >= h->size || left < 0)
        left = -1;
    if (right >= h->size || right < 0)
        right = -1;
 
    if (left != -1 && h->arr[left].distance < h->arr[index].distance)
        min = left;
    if (right != -1 && h->arr[right].distance < h->arr[min].distance)
        min = right;
 
    if (min != index) {
        struct heap_point temp = h->arr[min];
        h->arr[min] = h->arr[index];
        h->arr[index] = temp;
        heapify(h, min);
    }
}

struct point_distance_heap * build_heap(struct heap_point* data, int size, int capacity) {

    struct point_distance_heap* h = (struct point_distance_heap *) malloc(sizeof(struct point_distance_heap));
    if (h == NULL) {
        fprintf(stderr, "Error allocating memory for heap\n");
        return NULL;
    }

    h->capacity = capacity; 
    h->size = size;

    h->arr = data;

    int i = (h->size - 2) / 2;
    while (i >= 0) {
        heapify(h, i);
        i--;
    }
    return h;
}

void insert_helper(struct point_distance_heap* h, int index) {
 
    int parent = (index - 1) / 2;
 
    if (h->arr[parent].distance > h->arr[index].distance) {
        struct heap_point temp = h->arr[parent];
        h->arr[parent] = h->arr[index];
        h->arr[index] = temp;
        insert_helper(h, parent);
    }
}

struct heap_point heap_pop(struct point_distance_heap* h) {

    struct heap_point delete_item = {NULL, -1};
 
    if (h->size == 0) {
        printf("\nHeap id empty.");
        return delete_item;
    }
 
    delete_item = h->arr[0];
 
    h->arr[0] = h->arr[h->size - 1];
    h->size--;
 
    heapify(h, 0);
    return delete_item;
}
 
void heap_push(struct point_distance_heap* h, struct heap_point data) {
 
    // Checking if heap is full or not
    if (h->size < h->capacity) {
        // Inserting data into an array
        h->arr[h->size] = data;
        // Calling insertHelper function
        insert_helper(h, h->size);
        // Incrementing size of array
        h->size++;
    }
}
