#include "graph.h"

struct graph* new_graph () {
    struct graph *result = MALLOC_ONE_TYPE(struct graph);
    result->last_node = NULL;
    result->list_nodes = NULL;
    result->total_nodes = 0;
    result->total_edges = 0;
    result->has_point_data = false;

    return result;
}

void free_list_edges (struct node *n) {
    struct edge *e1 = n->list_edges;
    struct edge *e2 = n->list_edges->next;

    while (e1 != NULL) {
        free(e1);
        e1 = e2;
        if (e2 != NULL)
            e2 = e2->next;
    }

    n->list_edges = NULL;
}

void free_list_nodes (struct graph *g) {
    struct node *n1 = g->list_nodes;
    struct node *n2 = g->list_nodes->next;

    while (n1 != NULL) {
        if (n1->list_edges) {
            free_list_edges(n1);
        }

        free(n1);
        n1 = n2;
        if (n2 != NULL)
            n2 = n2->next;
    }
}

void free_graph (struct graph *g) {
    assert(g);

    if (g->list_nodes) {
        free_list_nodes(g);
    }

    free(g);
}

void insert_edge_graph (struct graph *g, const uint32_t id_1, const uint32_t id_2) {
    assert(g);

    struct node *n1, *n2;
    struct edge *edge;
    real_cpu norm;
    // Check if the edge is invalid
    if (id_1 == id_2) return;

    n1 = search_node(g,id_1);
    n2 = search_node(g,id_2);

    norm = calc_norm(n1->pos[0],n1->pos[1],n1->pos[2],n2->pos[0],n2->pos[1],n2->pos[2]);
    edge = new_edge(id_2,norm,n2);
    // First edge
    if (!n1->list_edges) {
        n1->list_edges = edge;
    } else { // Iterate over the list and insert to the last edge
        struct edge *e = n1->list_edges;
        while (e->next != NULL)
            e = e->next;
        e->next = edge;
    }
    // Increment the number of edges of origin Node
    n1->num_edges++;
    // Increment the total number of edges from the graph
    g->total_edges++;
}

void insert_node_graph (struct graph *g, const real_cpu pos[], const real_cpu sigma) {
    assert(g);

    struct node *tmp = g->list_nodes;
    struct node *node = new_node(g->total_nodes++,pos,sigma);
    // First node of the list
    if (!tmp) {
        g->list_nodes = node;
        g->last_node = node;
    } else { // Insert after the last node and update this pointer
        g->last_node->next = node;
        g->last_node = g->last_node->next;
    }
}

struct node* new_node (uint32_t id, const real_cpu pos[], const real_cpu sigma) {
    struct node *n = MALLOC_ONE_TYPE(struct node);
    n->id = id;
    memcpy(n->pos,pos,sizeof(real_cpu)*3);
    n->sigma = sigma;
    n->num_edges = 0;
    n->nmin_pmj = 10;           // Default values
    n->rpmj = 1000.0;           // Default values
    n->next = NULL;
    n->list_edges = NULL;

    return n;
}

struct edge* new_edge (uint32_t id, real_cpu w, struct node *dest) {
    struct edge *e = MALLOC_ONE_TYPE(struct edge);
    e->id = id;
    e->w = w;
    e->dest = dest;
    e->next = NULL;

    return e;
}

struct node* search_node (struct graph *g, const uint32_t id) {
    struct node *tmp = g->list_nodes;

    while (tmp != NULL) {
        if (tmp->id == id) {
            return tmp;
        }
        tmp = tmp->next;
    }

    fprintf(stderr,"[-] ERROR! Node %d was not found!\n",id);

    return NULL;
}

double* dijkstra (struct graph *g, const uint32_t src_id) {

    // Initialize the shortest distance array
    uint32_t num_nodes = g->total_nodes;
    double *dist = (double*)malloc(sizeof(double)*num_nodes);
    for (uint32_t i = 0; i < num_nodes; i++) dist[i] = __DBL_MAX__;
    dist[src_id] = 0.0;

    pqueue_t *pq;
    node_t   *ns;
    node_t   *n;

    // Initialize the priority queue
    ns = (struct node_t*)malloc(num_nodes * sizeof(node_t));
    pq = pqueue_init(num_nodes, cmp_pri, get_pri, set_pri, get_pos, set_pos);
    if (!(ns && pq)) {
        fprintf(stderr,"[graph] ERROR! Could not allocate priority queue!\n");
        exit(EXIT_FAILURE);
    }

    // Enqueue the source node
    ns[src_id].pri = 0.0;
    ns[src_id].val = src_id;
    pqueue_insert(pq, &ns[src_id]);

    while ((n = (node_t*)pqueue_pop(pq))) {

        double d = n->pri;
        uint32_t u = n->val;

        if (d > dist[u]) {
            continue;
        }

        struct node *u_node = search_node(g,u);
        struct edge *tmp = u_node->list_edges;

        while (tmp != NULL) {
            uint32_t v = tmp->id;
            double w = tmp->w;

            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;

                ns[v].pri = dist[v];
                ns[v].val = v;
                pqueue_insert(pq, &ns[v]);
            }

            tmp = tmp->next;
        }
    }

    pqueue_free(pq);
    free(ns);

    return dist;
}

void print_graph (struct graph *g) {

    struct node *n = g->list_nodes;

    while (n != NULL) {
        struct edge *e = n->list_edges;
        if (g->has_point_data) {
            fprintf(stdout,"|| %d (%.3lf,%.3lf,%.3lf) [%g] ||",n->id,n->pos[0],n->pos[1],n->pos[2],n->sigma);
        } else {
            fprintf(stdout,"|| %d (%.3lf,%.3lf,%.3lf) ||",n->id,n->pos[0],n->pos[1],n->pos[2]);
        }
        

        while (e != NULL) {
            fprintf(stdout," --> || %d %.3lf (%.3lf,%.3lf,%.3lf) ||",e->id,e->w,e->dest->pos[0],e->dest->pos[1],e->dest->pos[2]);
            e = e->next;
        }
        fprintf(stdout,"\n");

        n = n->next;
    }
    printf("Nodes = %u\n",g->total_nodes);
    printf("Edges = %u\n",g->total_edges);
    printf("Has point data = %d\n",g->has_point_data);
}

real_cpu calc_norm (const real_cpu x1, const real_cpu y1, const real_cpu z1,
                    const real_cpu x2, const real_cpu y2, const real_cpu z2) {
    return sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
}

bool is_terminal (const struct node *n) {
    return (n->num_edges == 1 && n->id != 0) ? true : false;
}
