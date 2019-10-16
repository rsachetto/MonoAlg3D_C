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
#include "../monodomain/constants.h"

struct node;
struct edge;

struct node
{
    uint32_t id;
    uint32_t num_edges;
    real_cpu x, y, z;

    struct edge *list_edges;
    struct node *next;
};

struct edge
{
    uint32_t id;
    real_cpu w;
    struct edge *next;
    struct node *dest;
};

struct graph
{
    struct node *list_nodes;
    struct node *last_node;
    uint32_t total_nodes;
    uint32_t total_edges;

    real_cpu dx;

    bool calc_retropropagation;
    real_cpu rpmj;
    real_cpu pmj_scale;
    uint32_t number_of_terminals;
};

struct node* new_node (uint32_t id, const real_cpu pos[]);
struct edge* new_edge (uint32_t id, real_cpu w, struct node *dest);
struct graph* new_graph ();

void insert_node_graph (struct graph *g, const real_cpu pos[]);
void insert_edge_graph (struct graph *g, const uint32_t id_1, const uint32_t id_2);
struct node* search_node (struct graph *g, const uint32_t id);

void print_graph (struct graph *g);
void free_graph (struct graph *g);
void free_list_nodes (struct graph *g);
void free_list_edges (struct node *n);

real_cpu calc_norm (const real_cpu x1, const real_cpu y1, const real_cpu z1,\
                  const real_cpu x2, const real_cpu y2, const real_cpu z2);

#endif //MONOALG3D_GRAPH_H
