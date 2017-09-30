//
// Created by sachetto on 29/09/17.
//

#include "cell.h"
#include <math.h>

void init_cell_with_cell_node(struct cell *data, bool init_ode) {

    data->type = CELL_NODE_TYPE;
    data->level   = 1;

    struct cell_node* cell_node = (struct cell_node*) malloc(sizeof(struct cell_node));
    init_cell_node(cell_node, init_ode);
    data->cell_node = cell_node;
    data->transition_node = NULL;

}

void init_cell_with_transition_node(struct cell *data) {

    data->type = TRANSITION_NODE_TYPE;
    data->level   = 1;

    struct transition_node* cell_node = (struct transition_node*) malloc(sizeof(struct transition_node));
    init_transition_node(cell_node);
    data->transition_node = cell_node;
    data->cell_node = NULL;

}

void init_cell_node(struct cell_node *cell_node, bool init_ode) {

    cell_node->center_x = 0.0;
    cell_node->center_y = 0.0;
    cell_node->center_z = 0.0;

    cell_node->active = true;

    cell_node->bunch_number = 0;

    cell_node->north       = NULL;
    cell_node->south       = NULL;
    cell_node->east        = NULL;
    cell_node->west        = NULL;
    cell_node->front       = NULL;
    cell_node->back        = NULL;

    cell_node->previous    = NULL;
    cell_node->next        = NULL;

    cell_node->grid_position        = 0;
    cell_node->gpu_sv_position      = 0;
    cell_node->hilbert_shape_number = 0;
    cell_node->face_length          = 1.0;
    cell_node->half_face_length     = 0.5;

    cell_node->v = 0;

    cell_node->north_flux = 0.0;
    cell_node->south_flux = 0.0;
    cell_node->east_flux  = 0.0;
    cell_node->west_flux  = 0.0;
    cell_node->front_flux = 0.0;
    cell_node->back_flux  = 0.0;

    cell_node->b = 0.0;

    cell_node->can_change = true;
    cell_node->fibrotic = false;
    cell_node->border_zone = false;
    cell_node->scar_type = 'n';


    //cell_node->firstElement = NULL; //TODO: @Incomplete

    //TODO: @Incomplete
    /*
    if(init_ode) {

        //firstElement = new Element; //TODO: @Incomplete

        od = (ode *) malloc(sizeof(ode));
        if (!od) {
            fprintf(stderr, "Error allocating memory for ode\n");
            exit(0);
        }

        //pthread_mutex_init(&updating, NULL);

    }
    else {
		od = NULL;
	}
     */

}

void free_cell_node(struct cell *cell) {

    //TODO: check cell type to free correctly

    if(cell->type == CELL_NODE_TYPE) {
        //TODO: @Incomplete
        /*
        if(od != NULL)
            free(od);

        if (firstElement) {
            Element *aux = firstElement;
            while(aux) {
                Element *temp = aux;
                aux = aux->next;
                delete temp;
            }
        }

        pthread_mutex_destroy(&updating);
         */
    }

    if(cell->cell_node)
        free(cell->cell_node);

    if(cell->transition_node)
        free(cell->transition_node);
};

void init_cell_node_ode(struct cell_node *cell_node) {


    //TODO: @Incomplete
    /*
     pthread_mutex_init(&updating, NULL);

	if (firstElement) {
		Element *aux = firstElement;
		while(aux) {
			Element *temp = aux;
			aux = aux->next;
			delete temp;
		}
	}

	firstElement = new Element;

    if(od != NULL)
        free(od);
    od = (ode *) malloc(sizeof(ode));
    if (!od) {
        std::cout << "Error allocating memory for ode\n" << std::endl;
        exit(0);
    }
     */
}

void lock_cell_node(struct cell_node *cell_node1) {

    //TODO: @Incomplete
    //pthread_mutex_lock(&(cell_node->updating));
}

void unlock_cell_node(struct cell_node *cell_node1) {

    //TODO: @Incomplete
    //pthread_mutex_unlock(&(cell_node->updating));
}

void init_transition_node(struct transition_node *transition_node) {

    transition_node->single_connector      = NULL;
    transition_node->quadruple_connector1  = NULL;
    transition_node->quadruple_connector2  = NULL;
    transition_node->quadruple_connector3  = NULL;
    transition_node->quadruple_connector4  = NULL;
    transition_node->direction             = 0;

}


void set_transition_node_data(struct cell *the_cell, char direction, struct cell *single_connector,
                              struct cell * quadruple_connector1, struct cell * quadruple_connector2,
                              struct cell * quadruple_connector3, struct cell * quadruple_connector4 ) {

    the_cell->transition_node->direction = direction;
    the_cell->transition_node->single_connector = single_connector;
    the_cell->transition_node->quadruple_connector1 = quadruple_connector1;
    the_cell->transition_node->quadruple_connector2 = quadruple_connector2;
    the_cell->transition_node->quadruple_connector3 = quadruple_connector3;
    the_cell->transition_node->quadruple_connector4 = quadruple_connector4;
}

void set_cell_node_data(struct cell *the_cell, float face_length, float half_face_length, uint64_t bunch_number,
                        struct cell *east, struct cell *north, struct cell *west, struct cell *south,
                        struct cell *front, struct cell *back,
                        struct cell *previous, struct cell *next,
                        uint64_t grid_position, uint8_t hilbert_shape_number,
                        float center_x, float center_y, float center_z)
{
    the_cell->cell_node->face_length = face_length;
    the_cell->cell_node->half_face_length = half_face_length;
    the_cell->cell_node->bunch_number = bunch_number;
    the_cell->cell_node->east = east;
    the_cell->cell_node->north = north;
    the_cell->cell_node->west = west;
    the_cell->cell_node->south = south;
    the_cell->cell_node->front = front;
    the_cell->cell_node->back = back;
    the_cell->cell_node->previous = previous;
    the_cell->cell_node->next = next;
    the_cell->cell_node->grid_position = grid_position;
    the_cell->cell_node->hilbert_shape_number = hilbert_shape_number;
    the_cell->cell_node->center_x = center_x;
    the_cell->cell_node->center_y = center_y;
    the_cell->cell_node->center_z = center_z;
}

void set_cell_flux( struct cell *the_cell, char direction ) {

    struct cell* neighbour_grid_cell;
    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    switch( direction ) {
        case 'n':
            neighbour_grid_cell = the_cell->cell_node->north;
            break;

        case 's':
            neighbour_grid_cell = the_cell->cell_node->south;
            break;

        case 'e':
            neighbour_grid_cell = the_cell->cell_node->east;
            break;

        case 'w':
            neighbour_grid_cell = the_cell->cell_node->west;
            break;

        case 'f':
            neighbour_grid_cell = the_cell->cell_node->front;
            break;

        case 'b':
            neighbour_grid_cell = the_cell->cell_node->back;
            break;
        default:
            fprintf(stderr, "Invalid cell direction %c! Exiting...", direction);
            exit(0);
    }

    float leastDistance = the_cell->cell_node->half_face_length;
    double localFlux;
    bool has_found;

    uint16_t neighbour_level =  neighbour_grid_cell->level;
    char neighbour_cell_type = neighbour_grid_cell->type;

    uint16_t the_cell_level = the_cell->level;

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    if (neighbour_level > the_cell_level ) {
        if((neighbour_cell_type == 'w') ) {
            white_neighbor_cell = neighbour_grid_cell->transition_node;
            has_found = false;
            while( !has_found ) {
                if( neighbour_cell_type == 'w' ) {
                    white_neighbor_cell = neighbour_grid_cell->transition_node;
                    if( white_neighbor_cell->single_connector == NULL ) {
                        has_found = true;
                    }
                    else {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_cell_type = neighbour_grid_cell->type;
                    }
                }
                else {
                    break;
                }
            }

        }
    }
        //Aqui, a célula vizinha tem um nivel de refinamento menor, entao eh mais simples.
    else {
        if(neighbour_level <= the_cell_level && (neighbour_cell_type == TRANSITION_NODE_TYPE) ) {
            white_neighbor_cell = neighbour_grid_cell->transition_node;
            has_found = false;
            while( !has_found ) {
                if( neighbour_cell_type == TRANSITION_NODE_TYPE ) {
                    white_neighbor_cell = neighbour_grid_cell->transition_node;
                    if( white_neighbor_cell->single_connector == NULL ) {
                        has_found = true;
                    }
                    else {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_cell_type = neighbour_grid_cell->type;
                    }
                }
                else {
                    break;
                }
            }
        }
    }

    neighbour_cell_type = neighbour_grid_cell->type;
    bool active = neighbour_grid_cell->cell_node->active;
    //Tratamos somente os pontos interiores da malha.
    if( ( neighbour_cell_type == CELL_NODE_TYPE ) && ( active == true ) )	{

        black_neighbor_cell = neighbour_grid_cell->cell_node;

        if ( black_neighbor_cell->half_face_length < leastDistance )
            leastDistance = black_neighbor_cell->half_face_length;

        localFlux = ( the_cell->cell_node->v - black_neighbor_cell->v ) / ( 2 * leastDistance );

        lock_cell_node(the_cell->cell_node);

        switch( direction ) {
            case 's':
                if ( localFlux > the_cell->cell_node->south_flux)
                    the_cell->cell_node->south_flux += localFlux;
                break;

            case 'n':
                if ( localFlux > the_cell->cell_node->north_flux)
                    the_cell->cell_node->north_flux += localFlux;
                break;

            case 'e':
                if ( localFlux > the_cell->cell_node->east_flux )
                    the_cell->cell_node->east_flux += localFlux;
                break;

            case 'w':
                if ( localFlux > the_cell->cell_node->west_flux)
                    the_cell->cell_node->west_flux += localFlux;
                break;

            case 'f':
                if ( localFlux > the_cell->cell_node->front_flux )
                    the_cell->cell_node->front_flux += localFlux;
                break;

            case 'b':
                if ( localFlux > the_cell->cell_node->back_flux )
                    the_cell->cell_node->back_flux += localFlux;
                break;
        }

        unlock_cell_node(the_cell->cell_node);

        lock_cell_node(black_neighbor_cell);

        switch( direction ) {
            case 's':
                if ( localFlux > black_neighbor_cell->north_flux )
                    black_neighbor_cell->north_flux += localFlux;
                break;

            case 'n':
                if ( localFlux > black_neighbor_cell->south_flux)
                    black_neighbor_cell->south_flux += localFlux;
                break;

            case 'e':
                if ( localFlux > black_neighbor_cell->west_flux)
                    black_neighbor_cell->west_flux += localFlux;
                break;

            case 'w':
                if ( localFlux > black_neighbor_cell->east_flux)
                    black_neighbor_cell->east_flux += localFlux;
                break;

            case 'f':
                if ( localFlux > black_neighbor_cell->back_flux)
                    black_neighbor_cell->back_flux += localFlux;
                break;

            case 'b':
                if ( localFlux > black_neighbor_cell->front_flux)
                    black_neighbor_cell->front_flux += localFlux;
                break;
        }

        unlock_cell_node(black_neighbor_cell);

    }
}

double get_cell_maximum_flux(struct cell_node* the_cell) {

    double maximumFlux = fabs(the_cell->east_flux);
    if( fabs(the_cell->north_flux) > maximumFlux )
        maximumFlux = fabs(the_cell->north_flux);

    if( fabs(the_cell->south_flux) > maximumFlux )
        maximumFlux = fabs(the_cell->south_flux);

    if( fabs(the_cell->west_flux) > maximumFlux )
        maximumFlux = fabs(the_cell->west_flux);

    if( fabs(the_cell->front_flux) > maximumFlux )
        maximumFlux = fabs(the_cell->front_flux);

    if( fabs(the_cell->back_flux) > maximumFlux )
        maximumFlux = fabs(the_cell->back_flux);

    return maximumFlux;
}

int getFreeSvPosition(short *gridToSV, int size) {
    for (int i = 0; i < size; ++i) {
        if(gridToSV[i] == -1) {
            gridToSV[i] = 0;
            return i;
        }
    }
    return 0;
}


void refine_cell( struct cell *cell, bool using_gpu, bool init_ode )  {
    if( cell == NULL ) {
        fprintf(stderr, "refine_cell(): Parameter cell is NULL. Exiting");
        exit(10);
    }

    size_t cell_size = sizeof(struct cell);

    struct cell_node* cell_node = cell->cell_node;

    struct cell *east_transition_cell,
            *north_transition_cell,
            *west_transition_cell,
            *south_transition_cell,
            *front_transition_cell,
            *back_transition_cell;

    struct cell *front_northeast_sub_cell,
            *front_northwest_sub_cell,
            *front_southwest_sub_cell,
            *front_southeast_sub_cell,
            *back_northeast_sub_cell,
            *back_northwest_sub_cell,
            *back_southwest_sub_cell,
            *back_southeast_sub_cell;

    int number_of_hilbert_shape;

    float cell_center_x       = cell->cell_node->center_x,
            cell_center_y     = cell->cell_node->center_y,
            cell_center_z     = cell->cell_node->center_z,
            cell_half_side    = cell_node->half_face_length,
            cell_quarter_side = cell_node->half_face_length / 2.0f;

    uint64_t old_bunch_number = cell_node->bunch_number;

    // Creation of the front northeast cell. This cell, which is to be refined,
    // becomes the frontNortheast cell of the new bunch.
    front_northeast_sub_cell                     = cell;
    front_northeast_sub_cell->level    = cell->level + (uint16_t )1;
    front_northeast_sub_cell->cell_node->face_length        = cell_half_side;
    front_northeast_sub_cell->cell_node->half_face_length   = cell_quarter_side;
    front_northeast_sub_cell->cell_node->center_x = cell_center_x + cell_quarter_side;
    front_northeast_sub_cell->cell_node->center_y = cell_center_y + cell_quarter_side;
    front_northeast_sub_cell->cell_node->center_z = cell_center_z + cell_quarter_side;
    front_northeast_sub_cell->cell_node->bunch_number       = old_bunch_number * 10 + 1;

    if(using_gpu) {
        //TODO: @Implement GPU code not implemented
        //refinedThisStep.push_back(front_northeast_sub_cell->gpuSVPosition);
    }

    // Creation of back Northeast node.
    back_northeast_sub_cell = (struct cell*) malloc(cell_size);
    init_cell_with_cell_node(back_northeast_sub_cell, init_ode);

    set_refined_cell_data(back_northeast_sub_cell,
          front_northeast_sub_cell,
          cell_half_side,
          cell_quarter_side,
          cell_center_x - cell_quarter_side,
          cell_center_y + cell_quarter_side,
          cell_center_z + cell_quarter_side,
          old_bunch_number * 10 + 2, using_gpu
    );


    // Creation of back Northwest node.
    back_northwest_sub_cell = (struct cell*) malloc(cell_size);
    init_cell_with_cell_node(back_northwest_sub_cell, init_ode);
    set_refined_cell_data(back_northwest_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x - cell_quarter_side,
                          cell_center_y - cell_quarter_side,
                          cell_center_z + cell_quarter_side,
                          old_bunch_number * 10 + 3, using_gpu
    );

    // Creation of front Northwest node.
    front_northwest_sub_cell = (struct cell*) malloc(cell_size);
    init_cell_with_cell_node(front_northwest_sub_cell, init_ode);
    set_refined_cell_data(front_northwest_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x + cell_quarter_side,
                          cell_center_y - cell_quarter_side,
                          cell_center_z + cell_quarter_side,
                          old_bunch_number * 10 + 4, using_gpu
    );


    // Creation of front Southwest node.
    front_southwest_sub_cell = (struct cell*) malloc(cell_size);
    init_cell_with_cell_node(front_southwest_sub_cell, init_ode);
    set_refined_cell_data(front_southwest_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x + cell_quarter_side,
                          cell_center_y - cell_quarter_side,
                          cell_center_z - cell_quarter_side,
                          old_bunch_number * 10 + 5, using_gpu
    );


    // Creation of back Southwest node.
    back_southwest_sub_cell = (struct cell*) malloc(cell_size);
    init_cell_with_cell_node(back_southwest_sub_cell, init_ode);
    set_refined_cell_data(back_southwest_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x - cell_quarter_side,
                          cell_center_y - cell_quarter_side,
                          cell_center_z - cell_quarter_side,
                          old_bunch_number * 10 + 6, using_gpu
    );



    // Creation of back Southeast node.
    back_southeast_sub_cell = (struct cell*) malloc(cell_size);
    init_cell_with_cell_node(back_southeast_sub_cell, init_ode);
    set_refined_cell_data(back_southeast_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x - cell_quarter_side,
                          cell_center_y + cell_quarter_side,
                          cell_center_z - cell_quarter_side,
                          old_bunch_number * 10 + 7, using_gpu
    );


    // Creation of front Southeast node.
    front_southeast_sub_cell = (struct cell*) malloc(cell_size);
    init_cell_with_cell_node(front_southeast_sub_cell, init_ode);
    set_refined_cell_data(front_southeast_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x + cell_quarter_side,
                          cell_center_y + cell_quarter_side,
                          cell_center_z - cell_quarter_side,
                          old_bunch_number * 10 + 8, using_gpu
    );

    // west transition node.
    west_transition_cell = (struct cell*)malloc(cell_size);
    init_cell_with_transition_node(west_transition_cell);
    set_refined_transition_node_data(west_transition_cell, front_northeast_sub_cell, 'w');


    // north transition node.
    north_transition_cell = (struct cell*)malloc(cell_size);
    init_cell_with_transition_node(north_transition_cell);
    set_refined_transition_node_data(north_transition_cell, front_northeast_sub_cell, 'n');

    // south transition node.
    south_transition_cell = (struct cell*)malloc(cell_size);
    init_cell_with_transition_node(south_transition_cell);
    set_refined_transition_node_data(south_transition_cell, front_northeast_sub_cell, 's');

    // east transition node.
    east_transition_cell = (struct cell*)malloc(cell_size);
    init_cell_with_transition_node(east_transition_cell);
    set_refined_transition_node_data(east_transition_cell, front_northeast_sub_cell, 'e');


    // front transition node.
    front_transition_cell = (struct cell*)malloc(cell_size);
    init_cell_with_transition_node(front_transition_cell);
    set_refined_transition_node_data(front_transition_cell, front_northeast_sub_cell, 'f');

    // back transitionition node.
    back_transition_cell = (struct cell*)malloc(cell_size);
    init_cell_with_transition_node(back_transition_cell);
    set_refined_transition_node_data(back_transition_cell, front_northeast_sub_cell, 'b');

    // Linking of new cell nodes and transition nodes.
    front_northeast_sub_cell->cell_node->north = north_transition_cell;
    front_northeast_sub_cell->cell_node->south = front_southeast_sub_cell;
    front_northeast_sub_cell->cell_node->east  = east_transition_cell;
    front_northeast_sub_cell->cell_node->west  = front_northwest_sub_cell;
    front_northeast_sub_cell->cell_node->front = front_transition_cell;
    front_northeast_sub_cell->cell_node->back  = back_northeast_sub_cell;

    back_northeast_sub_cell->cell_node->north = north_transition_cell;
    back_northeast_sub_cell->cell_node->south = back_southeast_sub_cell;
    back_northeast_sub_cell->cell_node->east  = east_transition_cell;
    back_northeast_sub_cell->cell_node->west  = back_northwest_sub_cell;
    back_northeast_sub_cell->cell_node->front = front_northeast_sub_cell;
    back_northeast_sub_cell->cell_node->back  = back_transition_cell;

    back_northwest_sub_cell->cell_node->north = north_transition_cell;
    back_northwest_sub_cell->cell_node->south = back_southwest_sub_cell;
    back_northwest_sub_cell->cell_node->east  = back_northeast_sub_cell;
    back_northwest_sub_cell->cell_node->west  = west_transition_cell;
    back_northwest_sub_cell->cell_node->front = front_northwest_sub_cell;
    back_northwest_sub_cell->cell_node->back  = back_transition_cell;

    front_northwest_sub_cell->cell_node->north = north_transition_cell;
    front_northwest_sub_cell->cell_node->south = front_southwest_sub_cell;
    front_northwest_sub_cell->cell_node->east  = front_northeast_sub_cell;
    front_northwest_sub_cell->cell_node->west  = west_transition_cell;
    front_northwest_sub_cell->cell_node->front = front_transition_cell;
    front_northwest_sub_cell->cell_node->back  = back_northwest_sub_cell;

    front_southwest_sub_cell->cell_node->north = front_northwest_sub_cell;
    front_southwest_sub_cell->cell_node->south = south_transition_cell;
    front_southwest_sub_cell->cell_node->east  = front_southeast_sub_cell;
    front_southwest_sub_cell->cell_node->west  = west_transition_cell;
    front_southwest_sub_cell->cell_node->front = front_transition_cell;
    front_southwest_sub_cell->cell_node->back  = back_southwest_sub_cell;

    back_southwest_sub_cell->cell_node->north = back_northwest_sub_cell;
    back_southwest_sub_cell->cell_node->south = south_transition_cell;
    back_southwest_sub_cell->cell_node->east  = back_southeast_sub_cell;
    back_southwest_sub_cell->cell_node->west  = west_transition_cell;
    back_southwest_sub_cell->cell_node->front = front_southwest_sub_cell;
    back_southwest_sub_cell->cell_node->back  = back_transition_cell;

    back_southeast_sub_cell->cell_node->north = back_northeast_sub_cell;
    back_southeast_sub_cell->cell_node->south = south_transition_cell;
    back_southeast_sub_cell->cell_node->east  = east_transition_cell;
    back_southeast_sub_cell->cell_node->west  = back_southwest_sub_cell;
    back_southeast_sub_cell->cell_node->front = front_southeast_sub_cell;
    back_southeast_sub_cell->cell_node->back  = back_transition_cell;

    front_southeast_sub_cell->cell_node->north = front_northeast_sub_cell;
    front_southeast_sub_cell->cell_node->south = south_transition_cell;
    front_southeast_sub_cell->cell_node->east  = east_transition_cell;
    front_southeast_sub_cell->cell_node->west  = front_southwest_sub_cell;
    front_southeast_sub_cell->cell_node->front = front_transition_cell;
    front_southeast_sub_cell->cell_node->back  = back_southeast_sub_cell;

    /* Connects the cell nodes with the transition nodes.  Quadruple  connectors
    1, 2, 3 and 4 are connected to neighbor cells  in  the  way  depicted below.

    This choice is made consistent with the one made at the function simplifyRef(),
    so that when two transition nodes of  same  level  connected  through  their
    single connectors are eliminated, the  subsequent  linking  of  corresponding
    quadruple connectors is correctly done.

                    front face           back face

                       ______              4______3
                     /|     /|            /|     /|
                   4/_|___3/ |           /_|____/ |
                   |  |___|__|          |  |___|__|
                   | /    | /           | /1   | /2
                   |/_____|/            |/_____|/
                   1      2
            ===========================================
                    west face           east face

                      3______               ______3
                     /|     /|            /|     /|
                   4/_|__ _/ |           /_|___4/ |
                   |  |___|__|          |  |___|__|
                   | /2   | /           | /    | /2
                   |/_____|/            |/_____|/
                   1                           1
            ===========================================
                    north face             south face

                      4______3              ______
                     /|     /|            /|     /|
                   1/_|___2/ |           /_|____/ |
                   |  |___|__|          |  |___|__|
                   | /    | /           | /4   | /3
                   |/_____|/            |/_____|/
                                        1      2
            ===========================================
     */

    // Front face.
    front_transition_cell->transition_node->quadruple_connector1 = front_southwest_sub_cell;
    front_transition_cell->transition_node->quadruple_connector2 = front_southeast_sub_cell;
    front_transition_cell->transition_node->quadruple_connector3 = front_northeast_sub_cell;
    front_transition_cell->transition_node->quadruple_connector4 = front_northwest_sub_cell;

    // Back face.
    back_transition_cell->transition_node->quadruple_connector1 = back_southwest_sub_cell;
    back_transition_cell->transition_node->quadruple_connector2 = back_southeast_sub_cell;
    back_transition_cell->transition_node->quadruple_connector3 = back_northeast_sub_cell;
    back_transition_cell->transition_node->quadruple_connector4 = back_northwest_sub_cell;

    // West face.
    west_transition_cell->transition_node->quadruple_connector1 = front_southwest_sub_cell;
    west_transition_cell->transition_node->quadruple_connector2 = back_southwest_sub_cell;
    west_transition_cell->transition_node->quadruple_connector3 = back_northwest_sub_cell;
    west_transition_cell->transition_node->quadruple_connector4 = front_northwest_sub_cell;

    // East face.
    east_transition_cell->transition_node->quadruple_connector1 = front_southeast_sub_cell;
    east_transition_cell->transition_node->quadruple_connector2 = back_southeast_sub_cell;
    east_transition_cell->transition_node->quadruple_connector3 = back_northeast_sub_cell;
    east_transition_cell->transition_node->quadruple_connector4 = front_northeast_sub_cell;

    // North face.
    north_transition_cell->transition_node->quadruple_connector1 = front_northwest_sub_cell;
    north_transition_cell->transition_node->quadruple_connector2 = front_northeast_sub_cell;
    north_transition_cell->transition_node->quadruple_connector3 = back_northeast_sub_cell;
    north_transition_cell->transition_node->quadruple_connector4 = back_northwest_sub_cell;

    // South face.
    south_transition_cell->transition_node->quadruple_connector1 = front_southwest_sub_cell;
    south_transition_cell->transition_node->quadruple_connector2 = front_southeast_sub_cell;
    south_transition_cell->transition_node->quadruple_connector3 = back_southeast_sub_cell;
    south_transition_cell->transition_node->quadruple_connector4 = back_southwest_sub_cell;



    // Linking bunch neighbor cells to the transition nodes just created.
    struct cell_node *neighbour_cell_node = NULL;
    struct transition_node *neighbour_transition_node = NULL;


    /*==========================================================================
                                WEST TRANSITION NODE

      Points the west neighboring cell to the transition node.

                                         t  t               B: Bunch
                                         | /                n: neighboring cell
                               n -- t -- B -- t             w: Transition node
                                       / |
                                      t  t

     *=========================================================================*/
    char node_type = west_transition_cell->transition_node->single_connector->type;
    if( node_type == 'b' ) {
        neighbour_cell_node = west_transition_cell->transition_node->single_connector->cell_node;
        neighbour_cell_node->east = west_transition_cell;
    }
    else if( node_type == TRANSITION_NODE_TYPE ) {
        neighbour_transition_node = west_transition_cell->transition_node->single_connector->transition_node;

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = west_transition_cell;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = west_transition_cell;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = west_transition_cell;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = west_transition_cell;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = west_transition_cell;
    }

    /*==========================================================================
                             NORTH TRANSITION NODE

      Points the north neighboring cell to the transition node.

                                      n
                                      |                  B: Bunch
                                      t  t               n: neighboring cell
                                      | /                t: Transition node
                                 t -- B -- t
                                    / |
                                   t  t

    ==========================================================================*/
    node_type = north_transition_cell->transition_node->single_connector->type;
    if( node_type == CELL_NODE_TYPE ) {
        neighbour_cell_node = north_transition_cell->transition_node->single_connector->cell_node;
        neighbour_cell_node->south = north_transition_cell;
    }
    else if( node_type == TRANSITION_NODE_TYPE )	{
        neighbour_transition_node = north_transition_cell->transition_node->single_connector->transition_node;

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = north_transition_cell;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = north_transition_cell;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = north_transition_cell;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = north_transition_cell;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = north_transition_cell;
    }

    /*==========================================================================
                                SOUTH TRANSITION NODE

      Points the south neighboring cell to the transition node.

                                      t  t               B: Bunch
                                      | /                n: neighboring cell
                                 t -- B -- t             w: Transition node
                                    / |
                                   t  t
                                      |
                                      n

    ==========================================================================*/
    node_type = south_transition_cell->transition_node->single_connector->type;
    if( node_type == CELL_NODE_TYPE ) {
        neighbour_cell_node = south_transition_cell->transition_node->single_connector->cell_node;
        neighbour_cell_node->north = south_transition_cell;
    }
    else if( node_type == TRANSITION_NODE_TYPE )	{
        neighbour_transition_node = south_transition_cell->transition_node->single_connector->transition_node;

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = south_transition_cell;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = south_transition_cell;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = south_transition_cell;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = south_transition_cell;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = south_transition_cell;
    }

    /*==========================================================================
                              EAST TRANSITION NODE

        Points the east neighboring cell to the transition node.

                                     t  t               B: Bunch
                                     | /                n: neighboring cell
                                t -- B -- t -- n        w: Transition node
                                   / |
                                  t  t

    ==========================================================================*/
    node_type = east_transition_cell->transition_node->single_connector->type;
    if( node_type == CELL_NODE_TYPE ) {
        neighbour_cell_node = east_transition_cell->transition_node->single_connector->cell_node;
        neighbour_cell_node->west = east_transition_cell;
    }
    else if( node_type == TRANSITION_NODE_TYPE ) {
        neighbour_transition_node = east_transition_cell->transition_node->single_connector->transition_node;

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = east_transition_cell;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = east_transition_cell;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = east_transition_cell;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = east_transition_cell;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = east_transition_cell;
    }

    /*==========================================================================
                                FRONT TRANSITION NODE

      Points the east neighboring cell to the transition node.

                                     t  t               B: Bunch
                                     | /                n: neighboring cell
                                t -- B -- t             t: TRANSITION node
                                   / |
                                  t  t
                                 /
                                n

    ==========================================================================*/
    node_type = front_transition_cell->transition_node->single_connector->type;
    if( node_type == CELL_NODE_TYPE ) {
        neighbour_cell_node = front_transition_cell->transition_node->single_connector->cell_node;
        neighbour_cell_node->back = front_transition_cell;
    }
    else if( node_type == TRANSITION_NODE_TYPE ) {
        neighbour_transition_node = front_transition_cell->transition_node->single_connector->transition_node;

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = front_transition_cell;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = front_transition_cell;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = front_transition_cell;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = front_transition_cell;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = front_transition_cell;
    }

    /*==========================================================================
                              BACK TRANSITION NODE

      Points the east neighboring cell to the transition node.

                                          n             B: Bunch
                                         /              n: neighboring cell
                                     t  t               t: Transition Node node
                                     | /
                                t -- B -- t
                                   / |
                                  t  t

    ==========================================================================*/
    node_type = back_transition_cell->transition_node->single_connector->type;
    if( node_type == CELL_NODE_TYPE ) {
        neighbour_cell_node = back_transition_cell->transition_node->single_connector->cell_node;
        neighbour_cell_node->front = back_transition_cell;
    }
    else if( node_type == TRANSITION_NODE_TYPE ) {
        neighbour_transition_node = back_transition_cell->transition_node->single_connector->transition_node;

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = back_transition_cell;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = back_transition_cell;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = back_transition_cell;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = back_transition_cell;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = back_transition_cell;
    }


    /*==========================================================================
                ORDERING OF CELL NODES THROUGH HILBERT'S CURVE
    ==========================================================================*/
    number_of_hilbert_shape = front_northeast_sub_cell->cell_node->hilbert_shape_number;

    if( number_of_hilbert_shape == 0 )	{
        /* Shape 0
                            _______
                           /      /      b: begin
                          /     b/       e: end
                          |  ______
                          | /     /
                          |/    e/
         */

        front_northeast_sub_cell->cell_node->hilbert_shape_number = 1;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 2;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 2;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 3;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 3;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 4;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 4;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 5;

        front_southeast_sub_cell->cell_node->next = front_northeast_sub_cell->cell_node->next;
        front_northeast_sub_cell->cell_node->next = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = front_southeast_sub_cell;

        front_southeast_sub_cell->cell_node->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = front_northeast_sub_cell;

        if( front_southeast_sub_cell->cell_node->next != 0 )
            front_southeast_sub_cell->cell_node->next->cell_node->previous = front_southeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 1 ) {
        /* Shape 1
                                       e
                               /|      |      b: begin
                             /  |   b  |      e: end
                             |  |___|__|
                             |      |
                             |______|
         */

        front_northeast_sub_cell->cell_node->hilbert_shape_number = 0;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 2;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 2;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 6;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 6;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 7;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 7;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 8;

        back_northeast_sub_cell->cell_node->next  = front_northeast_sub_cell->cell_node->next;
        front_northeast_sub_cell->cell_node->next = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = back_northeast_sub_cell;

        back_northeast_sub_cell->cell_node->previous  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = front_northeast_sub_cell;

        if( back_northeast_sub_cell->cell_node->next != 0 )
            back_northeast_sub_cell->cell_node->next->cell_node->previous = back_northeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 2 ) {
        /* Shape 2
                               /|     /|      b: begin
                             e/ |   b/ |      e: end
                                |      |
                               /      /
                              /______/
         */

        front_northeast_sub_cell->cell_node->hilbert_shape_number = 1;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 0;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 0;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 9;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 9;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 10;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 10;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 11;

        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell->cell_node->next;
        front_northeast_sub_cell->cell_node->next = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = front_northwest_sub_cell;

        front_northwest_sub_cell->cell_node->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = front_northeast_sub_cell;

        if( front_northwest_sub_cell->cell_node->next != 0 )
            front_northwest_sub_cell->cell_node->next->cell_node->previous = front_northwest_sub_cell;

    }

    else if( number_of_hilbert_shape == 3 ) {
        /* Shape 3
                               /b     /|      b: begin
                              /______/ |      e: end
                                       |
                               /e     /
                              /______/
         */

        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 11;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 7;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 7;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 0;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 0;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 9;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 9;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 6;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = back_northwest_sub_cell;

        back_southwest_sub_cell->cell_node->next  = front_northeast_sub_cell->cell_node->next;
        back_northwest_sub_cell->cell_node->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = back_southwest_sub_cell;

        back_northwest_sub_cell->cell_node->previous  = front_northeast_sub_cell->cell_node->previous;
        back_southwest_sub_cell->cell_node->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = back_northwest_sub_cell;

        if( back_southwest_sub_cell->cell_node->next != 0 )
            back_southwest_sub_cell->cell_node->next->cell_node->previous = back_southwest_sub_cell;

    }

    else if ( number_of_hilbert_shape == 4 ) {
        /* Shape 4
                                /|     /|      b: begin
                               /_|____/ |      e: end
                                 |      |
                                /      /
                              b/     e/
         */

        front_southwest_sub_cell->cell_node->hilbert_shape_number = 6;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 10;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 10;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 7;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 7;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 0;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 0;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 5;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = front_southwest_sub_cell;

        front_southeast_sub_cell->cell_node->next = front_northeast_sub_cell->cell_node->next;
        front_southwest_sub_cell->cell_node->next = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = front_southeast_sub_cell;

        front_southwest_sub_cell->cell_node->previous = front_northeast_sub_cell->cell_node->previous;
        front_southeast_sub_cell->cell_node->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = front_southwest_sub_cell;

        if ( front_southeast_sub_cell->cell_node->next != 0 )
            front_southeast_sub_cell->cell_node->next->cell_node->previous = front_southeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 5 ) {
        /* Shape 5
                                 ______
                                |      |      b: begin
                              __|___   |      e: end
                             |  |   |  |b
                             | /    |
                             |/     |e
         */

        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 8;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 9;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 9;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 11;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 11;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 4;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 4;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 0;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = back_southeast_sub_cell;

        front_southeast_sub_cell->cell_node->next = front_northeast_sub_cell->cell_node->next;
        back_southeast_sub_cell->cell_node->next  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = front_southeast_sub_cell;

        back_southeast_sub_cell->cell_node->previous  = front_northeast_sub_cell->cell_node->previous;
        front_southeast_sub_cell->cell_node->previous = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = back_southeast_sub_cell;

        if( front_southeast_sub_cell->cell_node->next != 0 )
            front_southeast_sub_cell->cell_node->next->cell_node->previous = front_southeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 6 ) {
        /* Shape 6
                                 ______
                                |      |      b: begin
                              __|___   |      e: end
                             | e|   |  |
                             |      | /
                            b|      |/
         */

        front_southwest_sub_cell->cell_node->hilbert_shape_number = 10;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 4;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 4;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 1;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 1;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 9;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 9;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 3;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = front_southwest_sub_cell;

        back_southwest_sub_cell->cell_node->next  = front_northeast_sub_cell->cell_node->next;
        front_southwest_sub_cell->cell_node->next = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = back_southwest_sub_cell;

        front_southwest_sub_cell->cell_node->previous = front_northeast_sub_cell->cell_node->previous;
        back_southwest_sub_cell->cell_node->previous  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = front_southwest_sub_cell;

        if( back_southwest_sub_cell->cell_node->next != 0 )
            back_southwest_sub_cell->cell_node->next->cell_node->previous = back_southwest_sub_cell;

    }

    else if( number_of_hilbert_shape == 7 ) {
        /* Shape 7
                                |b     |e     b: begin
                              __|___   |      e: end
                             |  |   |  |
                             | /    | /
                             |/     |/
         */

        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 3;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 11;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 11;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 4;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 4;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 1;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 1;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 8;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = back_northwest_sub_cell;

        back_northeast_sub_cell->cell_node->next  = front_northeast_sub_cell->cell_node->next;
        back_northwest_sub_cell->cell_node->next  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = back_northeast_sub_cell;

        back_northwest_sub_cell->cell_node->previous  = front_northeast_sub_cell->cell_node->previous;
        back_northeast_sub_cell->cell_node->previous  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = back_northwest_sub_cell;

        if( back_northeast_sub_cell->cell_node->next != 0 )
            back_northeast_sub_cell->cell_node->next->cell_node->previous = back_northeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 8 ) {
        /* Shape 8
                               /|     /e      b: begin
                              /_|____/        e: end
                                |
                               /      /b
                              /______/
         */

        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 5;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 9;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 9;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 10;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 10;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 7;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 7;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 1;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = back_southeast_sub_cell;

        back_northeast_sub_cell->cell_node->next  = front_northeast_sub_cell->cell_node->next;
        back_southeast_sub_cell->cell_node->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = back_northeast_sub_cell;

        back_southeast_sub_cell->cell_node->previous  = front_northeast_sub_cell->cell_node->previous;
        back_northeast_sub_cell->cell_node->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = back_southeast_sub_cell;

        if( back_northeast_sub_cell->cell_node->next != 0 )
            back_northeast_sub_cell->cell_node->next->cell_node->previous = back_northeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 9 ) {
        /* Shape 9
                                _______
                               /      /      b: begin
                              /      /       e: end
                              |      |
                              | /e   | /b
                              |/     |/
         */

        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 5;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 8;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 8;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 2;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 2;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 3;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 3;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 6;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = back_southeast_sub_cell;

        back_southwest_sub_cell->cell_node->next  = front_northeast_sub_cell->cell_node->next;
        back_southeast_sub_cell->cell_node->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->next = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = back_southwest_sub_cell;

        back_southeast_sub_cell->cell_node->previous  = front_northeast_sub_cell->cell_node->previous;
        back_southwest_sub_cell->cell_node->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->cell_node->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = back_southeast_sub_cell;

        if( back_southwest_sub_cell->cell_node->next != 0 )
            back_southwest_sub_cell->cell_node->next->cell_node->previous = back_southwest_sub_cell;

    }

    else if( number_of_hilbert_shape == 10 ) {
        /* Shape 10
                                 _______
                                /      /      b: begin
                              e/      /       e: end
                                  ____|__
                                 /    | /
                               b/     |/
         */

        front_southwest_sub_cell->cell_node->hilbert_shape_number = 6;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 4;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 4;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 8;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 8;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 2;
        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 2;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 11;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = front_southwest_sub_cell;

        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell->cell_node->next;
        front_southwest_sub_cell->cell_node->next = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->next  = front_northwest_sub_cell;

        front_southwest_sub_cell->cell_node->previous = front_northeast_sub_cell->cell_node->previous;
        front_northwest_sub_cell->cell_node->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = front_southwest_sub_cell;

        if( front_northwest_sub_cell->cell_node->next != 0 )
            front_northwest_sub_cell->cell_node->next->cell_node->previous = front_northwest_sub_cell;

    }

    else if( number_of_hilbert_shape == 11 ) {
        /* Shape 11
                                     b______
                                            |      b: begin
                                   e______  |      e: end
                                     _____|_|
                                    /     |
                                   /______|
         */

        back_northwest_sub_cell->cell_node->hilbert_shape_number  = 7;
        back_northeast_sub_cell->cell_node->hilbert_shape_number  = 3;
        back_southeast_sub_cell->cell_node->hilbert_shape_number  = 3;
        back_southwest_sub_cell->cell_node->hilbert_shape_number  = 5;
        front_southwest_sub_cell->cell_node->hilbert_shape_number = 5;
        front_southeast_sub_cell->cell_node->hilbert_shape_number = 10;
        front_northeast_sub_cell->cell_node->hilbert_shape_number = 10;
        front_northwest_sub_cell->cell_node->hilbert_shape_number = 2;

        if( front_northeast_sub_cell->cell_node->previous != 0 )
            front_northeast_sub_cell->cell_node->previous->cell_node->next = back_northwest_sub_cell;

        front_northwest_sub_cell->cell_node->next = front_northeast_sub_cell->cell_node->next;
        back_northwest_sub_cell->cell_node->next  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->next  = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->next  = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->next = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->next = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->next = front_northwest_sub_cell;

        back_northwest_sub_cell->cell_node->previous  = front_northeast_sub_cell->cell_node->previous;
        front_northwest_sub_cell->cell_node->previous = front_northeast_sub_cell;
        front_northeast_sub_cell->cell_node->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->cell_node->previous = front_southwest_sub_cell;
        front_southwest_sub_cell->cell_node->previous = back_southwest_sub_cell;
        back_southwest_sub_cell->cell_node->previous  = back_southeast_sub_cell;
        back_southeast_sub_cell->cell_node->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->cell_node->previous  = back_northwest_sub_cell;

        if( front_northwest_sub_cell->cell_node->next != 0 )
            front_northwest_sub_cell->cell_node->next->cell_node->previous = front_northwest_sub_cell;

    }

    // If necessary, simplifies the graph by eliminating adjacent transition nodes
    // of same level connected through their single connectors.
    simplify_refinement( east_transition_cell  );
    simplify_refinement( north_transition_cell );
    simplify_refinement( west_transition_cell  );
    simplify_refinement( south_transition_cell );
    simplify_refinement( front_transition_cell );
    simplify_refinement( back_transition_cell  );
}


/**
 * Simplifies data structure eliminating adjacent transition nodes of same level.
 *
 * @param cell_node Candidate transition node to be eliminated.
 *
 */
void simplify_refinement( struct cell *cell_node ) {

    if( cell_node == NULL ) {
       fprintf(stderr, "simplify_refinement: Parameter cell_node is NULL. Exiting!");
        exit(10);
    }

    // Pointers used to convert the Cell in a cell node or transition node.
    struct transition_node *neighbour_transition_node;
    struct cell_node *neighbour_cell_node;

    if( cell_node->transition_node->single_connector != 0 ) {

        // Both transition node and neighbor transition node must have the same
        // refinement level.
        char node_type = cell_node->type;
        uint16_t node_level = cell_node->level;

        uint16_t single_connector_level = cell_node->transition_node->single_connector->level;

        if( ( node_type == 'w') && (node_level == single_connector_level) ) {
            struct transition_node *neighbour_node = cell_node->transition_node->single_connector->transition_node;

            struct cell *cellNode[4];
            cellNode[0] = cell_node->transition_node->quadruple_connector1;
            cellNode[1] = cell_node->transition_node->quadruple_connector2;
            cellNode[2] = cell_node->transition_node->quadruple_connector3;
            cellNode[3] = cell_node->transition_node->quadruple_connector4;

            struct cell *neighborCell[4];
            neighborCell[0] = neighbour_node->quadruple_connector1;
            neighborCell[1] = neighbour_node->quadruple_connector2;
            neighborCell[2] = neighbour_node->quadruple_connector3;
            neighborCell[3] = neighbour_node->quadruple_connector4;

            char direction = cell_node->transition_node->direction;
            char type;

            for( int i = 0; i < 4; i++ ) {
                switch( direction ) {
                    case 'n': { cellNode[i]->cell_node->north = neighborCell[i]; break; }
                    case 's': { cellNode[i]->cell_node->south = neighborCell[i]; break; }
                    case 'e': { cellNode[i]->cell_node->east  = neighborCell[i]; break; }
                    case 'w': { cellNode[i]->cell_node->west  = neighborCell[i]; break; }
                    case 'f': { cellNode[i]->cell_node->front = neighborCell[i]; break; }
                    case 'b': { cellNode[i]->cell_node->back  = neighborCell[i]; break; }
                    default: break;
                }

                type = neighborCell[i]->type;
                switch( type ) {
                    case 'b': {
                        neighbour_cell_node = neighborCell[i]->cell_node;
                        switch( direction )	{
                            case 'n': { neighbour_cell_node->south = cellNode[i]; break; }
                            case 's': { neighbour_cell_node->north = cellNode[i]; break; }
                            case 'e': { neighbour_cell_node->west  = cellNode[i]; break; }
                            case 'w': { neighbour_cell_node->east  = cellNode[i]; break; }
                            case 'f': { neighbour_cell_node->back  = cellNode[i]; break; }
                            case 'b': { neighbour_cell_node->front = cellNode[i]; break; }
                            default: break;
                        }
                        break;
                    }

                    case 'w': {
                        neighbour_transition_node = neighborCell[i]->transition_node;
                        if( neighbour_node == neighbour_transition_node->single_connector->transition_node )
                            neighbour_transition_node->single_connector = cellNode[i];

                        else if( neighbour_node == neighbour_transition_node->quadruple_connector1->transition_node )
                            neighbour_transition_node->quadruple_connector1 = cellNode[i];

                        else if( neighbour_node == neighbour_transition_node->quadruple_connector2->transition_node)
                            neighbour_transition_node->quadruple_connector2 = cellNode[i];

                        else if( neighbour_node == neighbour_transition_node->quadruple_connector3->transition_node )
                            neighbour_transition_node->quadruple_connector3 = cellNode[i];

                        else if( neighbour_node == neighbour_transition_node->quadruple_connector4->transition_node )
                            neighbour_transition_node->quadruple_connector4 = cellNode[i];

                        break;
                    }

                    default: break;
                }
            }
            free_cell_node(cell_node);
            free(neighbour_node);
        }
    }
}

void set_refined_cell_data(struct cell* the_cell, struct cell* other_cell,
                           float face_length, float half_face_length,
                           float center_x, float center_y, float center_z,
                           uint64_t  bunch_number, bool using_gpu) {


    the_cell->level = other_cell->level;
    the_cell->cell_node->active = other_cell->cell_node->active;
    the_cell->cell_node->fibrotic = other_cell->cell_node->fibrotic;
    the_cell->cell_node->border_zone = other_cell->cell_node->border_zone;
    the_cell->cell_node->v = other_cell->cell_node->v;

    the_cell->cell_node->face_length = face_length;
    the_cell->cell_node->half_face_length = half_face_length;
    the_cell->cell_node->center_x = center_x;
    the_cell->cell_node->center_y = center_y;
    the_cell->cell_node->center_z = center_z;
    the_cell->cell_node->bunch_number = bunch_number;

    if (using_gpu) {
        //TODO: @Implement GPU related code
        /*
       the_cell->gpuSVPosition = freeSVPositions.back();
       freeSVPositions.pop_back();
       refinedThisStep.push_back(the_cell->gpuSVPosition);
        */
    }
    else {
        //TODO: @Incomplete: copy ODE stuff
        /*
        if(other_cell->od != NULL)
            memcpy (the_cell->od, other_cell->od, sizeof(ode));
        */
    }
}

void set_refined_transition_node_data(struct cell *the_node, struct cell* other_node, char direction) {


    the_node->transition_node->direction = direction;
    the_node->level                      = other_node->level;

    switch(direction) {
        case 'w': the_node->transition_node->single_connector   = other_node->cell_node->west; break;
        case 'n': the_node->transition_node->single_connector   = other_node->cell_node->north; break;
        case 's': the_node->transition_node->single_connector   = other_node->cell_node->south; break;
        case 'e': the_node->transition_node->single_connector   = other_node->cell_node->east; break;
        case 'f': the_node->transition_node->single_connector   = other_node->cell_node->front; break;
        case 'b': the_node->transition_node->single_connector   = other_node->cell_node->back; break;
        default:
            fprintf(stderr, "set_refined_transition_node_data() invalid direction %c Exiting!", direction);
            exit(10);
    }

}