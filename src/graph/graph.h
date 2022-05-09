//
// Created by bergolho on 20/07/18.
//

#ifndef MONOALG3D_GRAPH_H
#define MONOALG3D_GRAPH_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../common_types/common_types.h"

#include "pqueue.h"

struct node;
struct edge;

struct node {
    uint32_t id;
    uint32_t num_edges;
    uint32_t nmin_pmj;
    real_cpu sigma;
    real_cpu rpmj;
    real_cpu pos[3];

    struct edge *list_edges;
    struct node *next;
};

struct edge {
    uint32_t id;
    real_cpu w;
    struct edge *next;
    struct node *dest;
};

struct graph {
    struct node *list_nodes;
    struct node *last_node;
    uint32_t total_nodes;
    uint32_t total_edges;

    real_cpu dx;
    real_cpu rpmj;
    real_cpu pmj_scale;
    real_cpu asymm_ratio;
    uint32_t nmin_pmj;
    uint32_t nmax_pmj;
    uint32_t number_of_terminals;

    char *pmj_location_filename;

    bool has_pmj_location;
    bool has_point_data;
    bool calc_retropropagation;
};

struct node* new_node (uint32_t id, const real_cpu pos[], const real_cpu sigma);
struct edge* new_edge (uint32_t id, real_cpu w, struct node *dest);
struct graph* new_graph ();

void insert_node_graph (struct graph *g, const real_cpu pos[], const real_cpu sigma);
void insert_edge_graph (struct graph *g, const uint32_t id_1, const uint32_t id_2);
struct node* search_node (struct graph *g, const uint32_t id);

void print_graph (struct graph *g);
void free_graph (struct graph *g);
void free_list_nodes (struct graph *g);
void free_list_edges (struct node *n);

real_cpu calc_norm (const real_cpu x1, const real_cpu y1, const real_cpu z1,\
                  const real_cpu x2, const real_cpu y2, const real_cpu z2);

double* dijkstra (struct graph *g, const uint32_t src_id);

bool is_terminal (const struct node *n);

// --------------------------------------------------------------------------------
// PQUEUE library
typedef struct node_t {
    pqueue_pri_t pri;
    uint32_t     val;
    size_t       pos;
} node_t;


static int cmp_pri(pqueue_pri_t next, pqueue_pri_t curr) {
    return (next > curr);           // equivalent to std::greater<int>()
}


static pqueue_pri_t get_pri(void *a) {
    return ((node_t *) a)->pri;
}


static void set_pri(void *a, pqueue_pri_t pri) {
    ((node_t *) a)->pri = pri;
}


static size_t get_pos(void *a) {
    return ((node_t *) a)->pos;
}


static void set_pos(void *a, size_t pos) {
    ((node_t *) a)->pos = pos;
}


#endif //MONOALG3D_GRAPH_H
