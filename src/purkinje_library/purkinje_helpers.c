//
// Created by bergolho on 19/07/18.
//

#include "purkinje_helpers.h"

#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"
#include "../utils/utils.h"
#include "../string/sds.h"

#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef _MSC_VER
#include <process.h>
    #define getpid _getpid
#else
#include <unistd.h>
#endif

// Set a a custom Purkinje network from a file that stores its graph structure
void set_custom_purkinje_network (struct grid_purkinje *the_purkinje, const char *file_name, const real_cpu side_length, const real_cpu rpmj, const real_cpu pmj_scale, const bool calc_retro_propagation)
{

    struct graph *the_network = the_purkinje->the_network;

    set_purkinje_network_from_file(the_network,file_name,side_length,rpmj,pmj_scale,calc_retro_propagation);

    calculate_number_of_terminals(the_network);

    print_to_stdout_and_file("Number of Purkinje cells = %u\n",the_network->total_nodes);
    print_to_stdout_and_file("Number of Purkinje terminals = %u\n",the_network->number_of_terminals);

}

void set_purkinje_network_from_file (struct graph *the_purkinje_network, const char *file_name, const real_cpu side_length, const real_cpu rpmj, const real_cpu pmj_scale, const bool calc_retro_propagation)
{
    struct graph *skeleton_network = new_graph();

    //read_purkinje_network_from_file(file_name,&points,&branches,&N,&E);
    build_skeleton_purkinje(file_name,skeleton_network);

    build_mesh_purkinje(the_purkinje_network,skeleton_network,side_length,rpmj,pmj_scale,calc_retro_propagation);
    
    // Write the Purkinje to a VTK file for visualization purposes.
    write_purkinje_network_to_vtk(the_purkinje_network);
    //print_graph(the_purkinje_network);

    // Deallocate memory for the Skeleton mesh
    free_graph(skeleton_network);
    
}

// TODO: Find a way to build the purkinje network mesh directly without constructing the skeleton graph
void build_skeleton_purkinje (const char *filename, struct graph *skeleton_network)
{
    assert(skeleton_network);

    FILE *file = fopen(filename,"r");
    if (!file)
    {
        print_to_stdout_and_file("Error opening Purkinje mesh described in %s!!\n", filename);
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

void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *skeleton_network, const real_cpu side_length, const real_cpu rpmj, const real_cpu pmj_scale, const bool calc_retro_propagation)
{
    assert(the_purkinje_network);
    assert(skeleton_network);

    // TODO: Maybe avoid this conversion by building a skeleton mesh directly in micrometers
    // The side_length of a Purkinje volume is given in micrometers, so we convert to centimeters
    // um -> cm
    //the_purkinje_network->dx = side_length*UM_TO_CM;
    the_purkinje_network->dx = side_length;
    //the_purkinje_network->dx = side_length * 0.5;           // GAMBIARRA HARD !!!! Hermenegild
    the_purkinje_network->rpmj = rpmj;
    the_purkinje_network->pmj_scale = pmj_scale;
    the_purkinje_network->calc_retropropagation = calc_retro_propagation;

    uint32_t n = skeleton_network->total_nodes;
    // This map is needed to deal with bifurcations
    uint32_t *map_skeleton_to_mesh = (uint32_t*)calloc(n,sizeof(uint32_t));

    // Construct the first node
    struct node *tmp = skeleton_network->list_nodes;
    real_cpu pos[3]; pos[0] = tmp->x; pos[1] = tmp->y; pos[2] = tmp->z;
    insert_node_graph(the_purkinje_network,pos);
    
    // Make a Depth-First-Search to build the mesh of the Purkinje network
    depth_first_search(the_purkinje_network,tmp,0,map_skeleton_to_mesh);

    the_purkinje_network->dx = side_length;

    free(map_skeleton_to_mesh);

}

void depth_first_search (struct graph *the_purkinje_network, struct node *u, int level, uint32_t *map_skeleton_to_mesh)
{
    // TODO: Include the diameter of the Purkinje branch here ...

    struct edge *v = u->list_edges;
    while (v != NULL)
    {
        grow_segment(the_purkinje_network,u,v,map_skeleton_to_mesh);
        depth_first_search(the_purkinje_network,v->dest,level+1,map_skeleton_to_mesh);
        v = v->next;
    }

}

void grow_segment (struct graph *the_purkinje_network, struct node *u, struct edge *v, uint32_t *map_skeleton_to_mesh)
{
    real_cpu h = the_purkinje_network->dx;
    real_cpu d_ori[3], d[3];
    real_cpu segment_length = v->w;
    uint32_t n_points = segment_length / h;

    // Capture the index of the growing node on the mesh
    uint32_t id_source = map_skeleton_to_mesh[u->id];

    // Calculate a unitary direction vector of the segment
    calc_unitary_vector(d_ori,u,v->dest);

    // Copy the position of the source node
    d[0] = u->x;
    d[1] = u->y;
    d[2] = u->z;

    // DEBUG
    //print_to_stdout_and_file("Node %d will grow %d points\n",u->id,n_points);

    // Grow the number of points of size 'h' until reaches the size of the segment
    for (int k = 1; k <= n_points; k++)
    {
        real_cpu pos[3];
        pos[0] = d[0] + d_ori[0]*h*k;
        pos[1] = d[1] + d_ori[1]*h*k;
        pos[2] = d[2] + d_ori[2]*h*k;

        insert_node_graph(the_purkinje_network,pos);
        insert_edge_graph(the_purkinje_network,id_source,the_purkinje_network->total_nodes-1);
        insert_edge_graph(the_purkinje_network,the_purkinje_network->total_nodes-1,id_source);
        
        id_source = the_purkinje_network->total_nodes-1;
    }

    // Save the last inserted node index, in case this node generates offsprings
    map_skeleton_to_mesh[v->id] = id_source;
}

void calc_unitary_vector (real_cpu d_ori[], struct node *u, struct node *v)
{
    d_ori[0] = v->x - u->x;
    d_ori[1] = v->y - u->y;
    d_ori[2] = v->z - u->z;
    real_cpu norm = sqrt(d_ori[0]*d_ori[0] + d_ori[1]*d_ori[1] + d_ori[2]*d_ori[2]);
    for (int i = 0; i < 3; i++)
        d_ori[i] /= norm;
}

void read_purkinje_network_from_file (const char *filename, struct point **points, struct branch **branches, int *N, int *E)
{
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        print_to_stdout_and_file("Error opening Purkinje mesh described in %s!!\n", filename);
        exit (EXIT_FAILURE);
    }

    char str[100];

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    
    if (!fscanf(file,"%d",N))
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
    *points = (struct point*)malloc(sizeof(struct point)*(*N));
    for (int i = 0; i < (*N); i++)
    {
        real_cpu pos[3];
        if (!fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        } 
        (*points)[i].x = pos[0];
        (*points)[i].y = pos[1];
        (*points)[i].z = pos[2];
    }

    // Read edges
    int trash;
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"LINES") == 0) break;
    if (!fscanf(file,"%d %d",E,&trash))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    *branches = (struct branch*)malloc(sizeof(struct branch)*(*E));
    for (int i = 0; i < (*E); i++)
    {
        int e[2];
        if (!fscanf(file,"%d %d %d",&trash,&e[0],&e[1]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        }
        (*branches)[i].source = e[0];
        (*branches)[i].destination = e[1];
    }

    fclose(file);

}

// TODO: This function will be moved to the 'save_mesh_library'
void write_purkinje_network_to_vtk (struct graph *the_purkinje_network)
{
    assert(the_purkinje_network);

    struct node *n;
    struct edge *e;

    char *filename = "meshes/purkinje_mesh.vtk";
    print_to_stdout_and_file("Purkinje mesh file will be saved in :> %s\n",filename);

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

void calculate_number_of_terminals (struct graph *the_purkinje_network)
{
    assert(the_purkinje_network);

    uint32_t number_of_terminals = 0;

    struct node *n;
    struct edge *e;

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

// TODO: Some test for the network will be implemented here ...
int check_purkinje_input ()
{
    return 1;
}