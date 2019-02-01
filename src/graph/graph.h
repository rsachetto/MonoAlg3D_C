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

struct node;
struct edge;

struct node
{
    uint32_t id;
    uint32_t num_edges;
    double x, y, z;

    struct edge *list_edges;
    struct node *next;
};

struct edge
{
    uint32_t id;
    double w;
    struct edge *next;
    struct node *dest;
};

struct graph
{
    struct node *list_nodes;
    struct node *last_node;
    uint32_t total_nodes;
    uint32_t total_edges;

    double dx;
};

struct node* new_node (uint32_t id, const double pos[]);
struct edge* new_edge (uint32_t id, double w, struct node *dest);
struct graph* new_graph ();

void insert_node_graph (struct graph *g, const double pos[]);
void insert_edge_graph (struct graph *g, const uint32_t id_1, const uint32_t id_2);
struct node* search_node (struct graph *g, const uint32_t id);

void print_graph (struct graph *g);
void free_graph (struct graph *g);
void free_list_nodes (struct graph *g);
void free_list_edges (struct node *n);

double calc_norm (const double x1, const double y1, const double z1,\
                  const double x2, const double y2, const double z2);

#endif //MONOALG3D_GRAPH_H
