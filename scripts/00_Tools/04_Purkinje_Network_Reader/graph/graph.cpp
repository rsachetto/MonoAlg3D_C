#include "graph.h"

struct graph* new_graph ()
{
    struct graph *result = (struct graph*)malloc(sizeof(struct graph));
    result->last_node = NULL;
    result->list_nodes = NULL;
    result->total_nodes = 0;
    result->total_edges = 0;

    return result;
}

void free_list_edges (struct node *n)
{
    struct edge *e1 = n->list_edges;
    struct edge *e2 = n->list_edges->next;

    while (e1 != NULL)
    {
        e1->next = NULL;
        free(e1);
        e1 = e2;
        if (e2 != NULL)
            e2 = e2->next;
    }

    n->list_edges = NULL;
}

void free_list_nodes (struct graph *g)
{
    struct node *n1 = g->list_nodes;
    struct node *n2 = g->list_nodes->next;

    while (n1 != NULL)
    {
        if (n1->list_edges)
            free_list_edges(n1);
        
        n1->next = NULL;
        free(n1);
        n1 = n2;
        if (n2 != NULL)
            n2 = n2->next;
    }
}

void free_graph (struct graph *g)
{
    assert(g);

    if (g->list_nodes)
        free_list_nodes(g);
    
    free(g);
}

void insert_edge_graph (struct graph *g, const uint32_t id_1, const uint32_t id_2)
{
    assert(g);

    struct node *n1, *n2;
	struct edge *edge;
	double norm;
	// Check if the edge is invalid
	if (id_1 == id_2) return;

	n1 = search_node(g,id_1);
	n2 = search_node(g,id_2);
	
    norm = calc_norm(n1->x,n1->y,n1->z,n2->x,n2->y,n2->z);
    edge = new_edge(id_2,norm,n2);
    // First edge
    if (!n1->list_edges)
        n1->list_edges = edge;
    // Iterate over the list and insert to the last edge
    else
    {
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

void insert_node_graph (struct graph *g, const double pos[])
{
    assert(g);

    struct node *tmp = g->list_nodes;
    struct node *node = new_node(g->total_nodes++,pos);
    // First node of the list
    if (!tmp)
    {
        g->list_nodes = node;
        g->last_node = node;
    }
    // Insert after the last node and update this pointer
    else
    {
        g->last_node->next = node;
        g->last_node = g->last_node->next;
    }
}

void insert_node_graph (struct graph *g, const double pos[], const double value)
{
    assert(g);

    struct node *tmp = g->list_nodes;
    struct node *node = new_node(g->total_nodes++,pos,value);
    // First node of the list
    if (!tmp)
    {
        g->list_nodes = node;
        g->last_node = node;
    }
    // Insert after the last node and update this pointer
    else
    {
        g->last_node->next = node;
        g->last_node = g->last_node->next;
    }
}

struct node* new_node (uint32_t id, const double pos[])
{
    struct node *n = (struct node*)malloc(sizeof(struct node));
    n->id = id;
    n->x = pos[0];
    n->y = pos[1];
    n->z = pos[2];
    n->num_edges = 0;
    n->next = NULL;
    n->list_edges = NULL;
    
    return n;
}

struct node* new_node (uint32_t id, const double pos[], const double value)
{
    struct node *n = (struct node*)malloc(sizeof(struct node));
    n->id = id;
    n->x = pos[0];
    n->y = pos[1];
    n->z = pos[2];
    n->num_edges = 0;
    n->value = value;
    n->next = NULL;
    n->list_edges = NULL;
    
    return n;
}

struct edge* new_edge (uint32_t id, double w, struct node *dest)
{
    struct edge *e = (struct edge*)malloc(sizeof(struct edge));
	e->id = id;
	e->w = w;
	e->dest = dest;
	e->next = NULL;

	return e;
}

struct node* search_node (struct graph *g, const uint32_t id)
{
    struct node *tmp = g->list_nodes;
	while (tmp != NULL)
	{
		if (tmp->id == id)
			return tmp;
		tmp = tmp->next;
	}
    fprintf(stderr,"[-] ERROR! Node %d was not found!\n",id);
    
    return NULL;
}

void print_graph (struct graph *g)
{
    struct node *n = g->list_nodes;

    while (n != NULL)
    {
        struct edge *e = n->list_edges;
        fprintf(stdout,"|| %d (%.3lf,%.3lf,%.3lf) [%g] ||",n->id,n->x,n->y,n->z,n->value);

        while (e != NULL)
        {
            fprintf(stdout," --> || %d %.3lf (%.3lf,%.3lf,%.3lf) ||",e->id,e->w,e->dest->x,e->dest->y,e->dest->z);
            e = e->next;
        }
        fprintf(stdout,"\n");

        n = n->next;
    } 
    printf("Nodes = %u\n",g->total_nodes);
    printf("Edges = %u\n",g->total_edges);


}

double calc_norm (const double x1, const double y1, const double z1,\
                  const double x2, const double y2, const double z2)
{
    return sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
}

void build_graph_from_purkinje_network (struct graph *the_graph, struct purkinje_network *the_purkinje_network)
{
    uint32_t np = the_purkinje_network->num_points;
    uint32_t nl = the_purkinje_network->num_lines;

    struct point *the_points = the_purkinje_network->points;
    struct line *the_lines = the_purkinje_network->lines;
    double *the_data = the_purkinje_network->point_data;

    for (uint32_t i = 0; i < np; i++)
    {
        double pos[3];
        pos[0] = the_points[i].x;
        pos[1] = the_points[i].y;
        pos[2] = the_points[i].z;

        double value;
        value = the_data[i];

        insert_node_graph(the_graph,pos,value);
    }
        
    for (uint32_t i = 0; i < nl; i++)
    {
        uint32_t link[2];
        link[0] = the_lines[i].src;
        link[1] = the_lines[i].dest;
        
        insert_edge_graph(the_graph,link[0],link[1]);
    }

    //print_graph(the_graph);
}

struct terminal* calculate_terminals (struct graph *the_graph, uint32_t &num_terminals)
{
    num_terminals = count_number_of_terminals(the_graph);

    struct terminal *the_terminals = (struct terminal*)malloc(sizeof(struct terminal)*num_terminals);

    uint32_t i = 0;
    struct node *tmp = the_graph->list_nodes;
    while (tmp != NULL)
    {
        if (is_terminal(tmp))
        {
            the_terminals[i].id = tmp->id;
            the_terminals[i].x = tmp->x;
            the_terminals[i].y = tmp->y;
            the_terminals[i].z = tmp->z;
            the_terminals[i].value = tmp->value;

            i++;
        }

        tmp = tmp->next;
    }

    return the_terminals;
}

uint32_t count_number_of_terminals (struct graph *the_graph)
{
    uint32_t ret = 0;
    struct node *tmp = the_graph->list_nodes;

    while (tmp != NULL)
    {
        if (is_terminal(tmp))
            ret++;

        tmp = tmp->next;
    }

    return ret;
}

bool is_terminal (struct node *the_node)
{
    if (the_node->num_edges == 1 && the_node->id != 0)
        return true;
    else
        return false;
}

void write_configuration_file (struct terminal *the_terminals, const uint32_t num_terminals)
{
    char *filename = "outputs/lucas_cheat_purkinje_coupled.ini";
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads = %u\n",NUM_THREADS);
    fprintf(file,"dt_pde = %g\n",DT_PDE);
    fprintf(file,"simulation_time = %g\n",SIMULATION_TIME);
    if (ABORT_ON_NO_ACTIVITY)
        fprintf(file,"abort_on_no_activity = true\n");
    else
        fprintf(file,"abort_on_no_activity = false\n");
    if (USE_ADAPTIVITY)
        fprintf(file,"use_adaptivity = true\n");
    else
        fprintf(file,"use_adaptivity = false\n");
    fprintf(file,"\n");
    
    fprintf(file,"[save_result]\n");
    fprintf(file,"print_rate = %u\n",PRINT_RATE);
    fprintf(file,"output_dir = %s\n",OUTPUT_DIR);
    fprintf(file,"main_function = %s\n",SAVE_MESH_FUNCTION);
    if (SAVE_PVD)
        fprintf(file,"save_pvd = true\n");
    else
        fprintf(file,"save_pvd = false\n");
    fprintf(file,"file_prefix = %s\n",FILE_PREFIX);
    fprintf(file,"\n");

    fprintf(file,"[update_monodomain]\n");
    fprintf(file,"main_function = %s\n",UPDATE_MONODOMAIN_FUNCTION);
    fprintf(file,"\n");

    fprintf(file,"[assembly_matrix]\n");
    fprintf(file,"init_function=set_initial_conditions_fvm\n");
    fprintf(file,"sigma_x = %g\n",SIGMA_X);
    fprintf(file,"sigma_y = %g\n",SIGMA_Y);
    fprintf(file,"sigma_z = %g\n",SIGMA_Z);
    fprintf(file,"library_file = %s\n",ASSEMBLY_MATRIX_LIBRARY);
    fprintf(file,"main_function = %s\n",ASSEMBLY_MATRIX_FUNCTION);
    fprintf(file,"\n");

    fprintf(file,"[linear_system_solver]\n");
    fprintf(file,"tolerance = %g\n",TOLERANCE);
    if (USE_PRECONDITIONER)
        fprintf(file,"use_preconditioner = true\n");
    else
        fprintf(file,"use_preconditioner = false\n");
    fprintf(file,"max_iterations = %u\n",MAX_ITERATIONS);
    fprintf(file,"library_file = %s\n",LINEAR_SYSTEM_LIBRARY);
    fprintf(file,"main_function = %s\n",LINEAR_SYSTEM_FUNCTION);
    fprintf(file,"\n");

    fprintf(file,"[alg]\n");
    fprintf(file,"refinement_bound = %g\n",REFINEMENT_BOUND);
    fprintf(file,"derefinement_bound = %g\n",DEREFINEMENT_BOUND);
    fprintf(file,"refine_each = %u\n",REFINE_EACH);
    fprintf(file,"derefine_each = %u\n",DEREFINE_EACH);
    fprintf(file,"\n");

    fprintf(file,"[domain]\n");
    fprintf(file,"name = %s\n",DOMAIN_NAME);
    fprintf(file,"num_layers = %u\n",NUM_LAYER);
    fprintf(file,"start_dx = %g\n",START_DX);
    fprintf(file,"start_dy = %g\n",START_DY);
    fprintf(file,"start_dz = %g\n",START_DZ);
    fprintf(file,"side_length = %g\n",SIDE_LENGTH);
    fprintf(file,"main_function = %s\n",DOMAIN_FUNCTION);
    fprintf(file,"\n");

    fprintf(file,"[ode_solver]\n");
    fprintf(file,"dt_ode = %g\n",DT_ODE);
    if (USE_GPU)
        fprintf(file,"use_gpu = true\n");
    else
        fprintf(file,"use_gpu = false\n");
    fprintf(file,"gpu_id = %u\n",GPU_ID);
    fprintf(file,"library_file = %s\n",ODE_LIBRARY);
    fprintf(file,"\n");

    for (uint32_t i = 0; i < num_terminals; i++)
    {
        fprintf(file,"[stim_s%u]\n",i+1);
        fprintf(file,"start = %g\n",the_terminals[i].value);
        fprintf(file,"duration = %g\n",STIM_DURATION);
        fprintf(file,"current = %g\n",STIM_CURRENT);
        fprintf(file,"period = %g\n",STIM_PERIOD);
        fprintf(file,"center_x = %g\n",the_terminals[i].x);
        fprintf(file,"center_y = %g\n",the_terminals[i].y);
        fprintf(file,"center_z = %g\n",the_terminals[i].z);
        fprintf(file,"radius = %g\n",PMJ_RADIUS);
        fprintf(file,"main_function = %s\n",STIM_FUNCTION);
        fprintf(file,"\n");
    }

    fclose(file);

    printf("[!] Configuration file saved at :> %s\n",filename);
}