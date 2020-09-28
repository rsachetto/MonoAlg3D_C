//
// Created by bergolho on 19/07/18.
//

#include "purkinje_helpers.h"

#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"
#include "../utils/utils.h"
#include "../3dparty/sds/sds.h"

#include <math.h>
#include <string.h>

// Set a a custom Purkinje network from a file that stores its graph structure
void set_custom_purkinje_network (struct grid_purkinje *the_purkinje, const char *file_name, const real_cpu dx) {

    struct graph *the_network = the_purkinje->network;

    set_purkinje_network_from_file(the_network,file_name,dx);

    calculate_number_of_terminals(the_network);

    log_to_stdout_and_file("Number of Purkinje cells = %u\n",the_network->total_nodes);
    log_to_stdout_and_file("Number of Purkinje terminals = %u\n",the_network->number_of_terminals);

}

void set_purkinje_network_from_file (struct graph *the_purkinje_network, const char *file_name, const real_cpu dx) {

    struct graph *skeleton_network = new_graph();

    //read_purkinje_network_from_file(file_name,&points,&branches,&N,&E);
    build_skeleton_purkinje(file_name,skeleton_network);

    build_mesh_purkinje(the_purkinje_network,skeleton_network,dx);
    
    // Write the Purkinje to a VTK file for visualization purposes.
    write_purkinje_network_to_vtk(the_purkinje_network);
    //print_graph(the_purkinje_network);

    // Deallocate memory for the Skeleton mesh
    free_graph(skeleton_network);
    
}

void build_skeleton_purkinje (const char *filename, struct graph *skeleton_network)
{
    assert(skeleton_network);

    FILE *file = fopen(filename,"r");
    if (!file)
    {
        log_to_stdout_and_file("Error opening Purkinje mesh described in %s!!\n", filename);
        exit (EXIT_FAILURE);
    }

    int N;
    char str[100];

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    
    if (!fscanf(file,"%d",&N))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    if (!fscanf(file,"%s",str))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }

    // Read points
    for (int i = 0; i < N; i++)
    {
        real_cpu pos[3];
        if (!fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        } 
        insert_node_graph(skeleton_network,pos);
    }

    // Read edges
    int trash, E;
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"LINES") == 0) break;
    if (!fscanf(file,"%d %d",&E,&trash))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < E; i++)
    {
        int e[2];
        if (!fscanf(file,"%d %d %d",&trash,&e[0],&e[1]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        }
        insert_edge_graph(skeleton_network,e[0],e[1]);
    }

    fclose(file);
}

void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *skeleton_network, const real_cpu dx) {
    
    assert(the_purkinje_network);
    assert(skeleton_network);

    the_purkinje_network->dx = dx;

    uint32_t n = skeleton_network->total_nodes;
    // This map is needed to deal with bifurcations
    uint32_t *map_skeleton_to_mesh = (uint32_t*)calloc(n,sizeof(uint32_t));

    // Construct the first node
    struct node *tmp = skeleton_network->list_nodes;
    real_cpu pos[3], sigma; 
    pos[0] = tmp->x; 
    pos[1] = tmp->y; 
    pos[2] = tmp->z;
    insert_node_graph(the_purkinje_network,pos);
    
    // Initialize the DFS mask for the visited nodes
    bool *dfs_visited = (bool*)malloc(sizeof(bool)*n);
    for (uint32_t i = 0; i < n; i++) dfs_visited[i] = false;

    // Make a Depth-First-Search to build the mesh of the Purkinje network
    depth_first_search(the_purkinje_network,tmp,0,map_skeleton_to_mesh,dfs_visited);
    
    free(dfs_visited);
    free(map_skeleton_to_mesh);

}

void depth_first_search (struct graph *the_purkinje_network, struct node *u, int level, uint32_t *map_skeleton_to_mesh, bool *dfs_visited) {
    
    dfs_visited[u->id] = true;

    struct edge *v = u->list_edges;
    while (v != NULL)
    {
        if (dfs_visited[v->id] == false)
        {
            grow_segment(the_purkinje_network,u,v,map_skeleton_to_mesh);
            depth_first_search(the_purkinje_network,v->dest,level+1,map_skeleton_to_mesh,dfs_visited);
        }
        v = v->next;
    }

}

void grow_segment (struct graph *the_purkinje_network, struct node *u, struct edge *v, uint32_t *map_skeleton_to_mesh)
{
    real_cpu dx = the_purkinje_network->dx;
    real_cpu d_ori[3], d[3];
    real_cpu segment_length = v->w;
    real_cpu remainder_points = fmod(segment_length,dx);
    uint32_t n_points = segment_length / dx;

    // Capture the index of the growing node on the mesh
    uint32_t id_source = map_skeleton_to_mesh[u->id];

    // Calculate a unitary direction vector of the segment
    calc_unitary_vector(d_ori,u,v->dest);

    // Copy the position of the source node
    d[0] = u->x;
    d[1] = u->y;
    d[2] = u->z;

    // DEBUG
    //log_to_stdout_and_file("Node %d will grow %d points\n",u->id,n_points);

    // Grow the number of points of size 'h' until reaches the size of the segment
    for (int k = 1; k <= n_points; k++) {

        real_cpu pos[3];
        pos[0] = d[0] + d_ori[0]*dx*k;
        pos[1] = d[1] + d_ori[1]*dx*k;
        pos[2] = d[2] + d_ori[2]*dx*k;

        insert_node_graph(the_purkinje_network,pos);
        insert_edge_graph(the_purkinje_network,id_source,the_purkinje_network->total_nodes-1);
        insert_edge_graph(the_purkinje_network,the_purkinje_network->total_nodes-1,id_source);
        
        id_source = the_purkinje_network->total_nodes-1;
    }
    // Grow any remainder point with a size 'dx'
    if (remainder_points > 0.0) {

        real_cpu pos[3];
        pos[0] = d[0] + d_ori[0]*dx*n_points + d_ori[0]*remainder_points;
        pos[1] = d[1] + d_ori[1]*dx*n_points + d_ori[1]*remainder_points;
        pos[2] = d[2] + d_ori[2]*dx*n_points + d_ori[2]*remainder_points;

        insert_node_graph(the_purkinje_network,pos);
        insert_edge_graph(the_purkinje_network,id_source,the_purkinje_network->total_nodes-1);
        insert_edge_graph(the_purkinje_network,the_purkinje_network->total_nodes-1,id_source);
        
        id_source = the_purkinje_network->total_nodes-1;
    }

    // Save the last inserted node index, in case this node generates offsprings
    map_skeleton_to_mesh[v->id] = id_source;
}

void calc_unitary_vector (real_cpu d_ori[], struct node *u, struct node *v) {

    d_ori[0] = v->x - u->x;
    d_ori[1] = v->y - u->y;
    d_ori[2] = v->z - u->z;
    real_cpu norm = sqrt(d_ori[0]*d_ori[0] + d_ori[1]*d_ori[1] + d_ori[2]*d_ori[2]);
    for (int i = 0; i < 3; i++)
        d_ori[i] /= norm;
}


void write_purkinje_network_to_vtk (struct graph *the_purkinje_network) {

    assert(the_purkinje_network);

    struct node *n;
    struct edge *e;

    char *filename = "meshes/purkinje_mesh.vtk";
    log_to_stdout_and_file("Purkinje mesh file will be saved in :> %s\n",filename);

    FILE *file = fopen(filename,"w+");
    fprintf(file,"# vtk DataFile Version 3.0\n");
    fprintf(file,"Purkinje\nASCII\nDATASET POLYDATA\n");
    fprintf(file,"POINTS %d float\n",the_purkinje_network->total_nodes);

    n = the_purkinje_network->list_nodes;
    while (n != NULL)
    {
        fprintf(file,"%.8lf %.8lf %.8lf\n",n->x,n->y,n->z);
        n = n->next;
    }

    n = the_purkinje_network->list_nodes;
    fprintf(file,"LINES %d %d\n",the_purkinje_network->total_edges,the_purkinje_network->total_edges*3);
    while (n != NULL)
    {
        e = n->list_edges;
        while (e != NULL)
        {
            fprintf(file,"2 %d %d\n",n->id,e->id);
            e = e->next;
        }
        n = n->next;
    }

    fclose(file);
}

void calculate_number_of_terminals (struct graph *the_purkinje_network) {

    assert(the_purkinje_network);

    uint32_t number_of_terminals = 0;

    struct node *n;
    //struct edge *e;

    n = the_purkinje_network->list_nodes;
    while (n != NULL)
    {
        if (is_terminal(n))
            number_of_terminals++;

        n = n->next;
    }

    the_purkinje_network->number_of_terminals = number_of_terminals;
}

bool is_terminal (const struct node *n)
{
    if (n->num_edges == 1 && n->id != 0)
        return true;
    else
        return false;
}

// Check if there are duplicates points inside the Purkinje network
int check_purkinje_mesh_for_errors (struct graph *the_purkinje_network) {
    
    int no_duplicates = 1;

    // Check duplicates
    struct node *tmp = the_purkinje_network->list_nodes;
    while (tmp != NULL)
    {
        struct node *tmp2 = the_purkinje_network->list_nodes;
        while (tmp2 != NULL)
        {
            if (tmp->x == tmp2->x && tmp->y == tmp2->y && tmp->z == tmp2->z && tmp->id != tmp2->id)
            {
                printf("[purkinje] Duplicates are indexes: %u and %u --> (%g,%g,%g) x (%g,%g,%g)\n",tmp->id,tmp2->id,tmp->x,tmp->y,tmp->z,tmp2->x,tmp2->y,tmp2->z);
                printf("\t|| %u has %u edges || %u has %u edges ||\n",tmp->id,tmp->num_edges,tmp2->id,tmp2->num_edges);
                no_duplicates = 0;
            }
            tmp2 = tmp2->next;
        }
        tmp = tmp->next;
    }
    return no_duplicates;
}