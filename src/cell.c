//
// Created by sachetto on 29/09/17.
//

#include "cell.h"



void init_basic_cell_data(struct basic_cell_data *data, uint8_t level, float center_x, float center_y, float center_z) {

    data->type = 'x';
    data->level   = level;
    data->center_x = center_x;
    data->center_y = center_y;
    data->center_z = center_z;

}

void init_basic_cell_data_with_default_values(struct basic_cell_data *data, char type) {

    data->type = type;
    data->level   = 1;
    data->center_x = 0.0;
    data->center_y = 0.0;
    data->center_z = 0.0;

}

void init_cell_node(struct cell_node *cell_node, bool init_ode) {

    init_basic_cell_data_with_default_values(&(cell_node->cell_data), CELL_NODE_TYPE);

    cell_node->active = 1;

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

void free_cell_node(struct cell_node *cell_node) {

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

    init_basic_cell_data_with_default_values(&(transition_node->cell_data), TRANSITION_NODE_TYPE);

    transition_node->single_connector      = NULL;
    transition_node->quadruple_connector1  = NULL;
    transition_node->quadruple_connector2  = NULL;
    transition_node->quadruple_connector3  = NULL;
    transition_node->quadruple_connector4  = NULL;
    transition_node->direction             = 0;

}


void set_transition_node_data(struct transition_node *the_transtion_node, char direction, float center_x,
                              float center_y, float center_z, void *single_connector,
                              void * quadruple_connector1, void * quadruple_connector2,
                              void * quadruple_connector3, void * quadruple_connector4 ) {

    the_transtion_node->direction = direction;
    the_transtion_node->cell_data.center_x = center_x;
    the_transtion_node->cell_data.center_y = center_y;
    the_transtion_node->cell_data.center_z = center_z;
    the_transtion_node->single_connector = single_connector;

    the_transtion_node->quadruple_connector1 = quadruple_connector1;

    the_transtion_node->quadruple_connector2 = quadruple_connector2;

    the_transtion_node->quadruple_connector3 = quadruple_connector3;

    the_transtion_node->quadruple_connector4 = quadruple_connector4;
}

void set_cell_node_data(struct cell_node *the_cell, float face_length, float half_face_length, uint64_t bunch_number,
                        void *east, void *north, void *west, void *south, void *front, void *back,
                        void *previous, void *next,
                        uint64_t grid_position, uint8_t hilbert_shape_number,
                        float center_x, float center_y, float center_z)
{
    the_cell->face_length = face_length;
    the_cell->half_face_length = half_face_length;
    the_cell->bunch_number = bunch_number;
    the_cell->east = east;
    the_cell->north = north;
    the_cell->west = west;
    the_cell->south = south;
    the_cell->front = front;
    the_cell->back = back;
    the_cell->previous = previous;
    the_cell->next = next;
    the_cell->grid_position = grid_position;
    the_cell->hilbert_shape_number = hilbert_shape_number;
    the_cell->cell_data.center_x = center_x;
    the_cell->cell_data.center_y = center_y;
    the_cell->cell_data.center_z = center_z;
}

//void set_cell_flux( struct cell_node *the_cell, char direction ) {
//
//    union cell neighbor_grid_cell;
//    struct transition_node *white_neighbor_cell;
//    struct cell_node *black_neighbor_cell;
//
//    switch( direction ) {
//        case 'n':
//            neighbor_grid_cell = the_cell->north;
//            break;
//
//        case 's':
//            neighbor_grid_cell = the_cell->south;
//            break;
//
//        case 'e':
//            neighbor_grid_cell = the_cell->east;
//            break;
//
//        case 'w':
//            neighbor_grid_cell = the_cell->west;
//            break;
//
//        case 'f':
//            neighbor_grid_cell = the_cell->front;
//            break;
//
//        case 'b':
//            neighbor_grid_cell = the_cell->back;
//            break;
//
//    }
//
//    float leastDistance = the_cell->half_face_length;
//    double localFlux;
//    bool hasFound;
//
//
//    /* When neighbor_grid_cell is a transition node, looks for the next neighbor
//     * cell which is a cell node. */
//    //Acha uma célula real que está no caixo enviado como vizinho
//    if (neighbor_grid_cell->level > level ) {
//        if((neighbor_grid_cell->type == 'w') ) {
//            white_neighbor_cell = static_cast<TransitionNode*>(neighbor_grid_cell);
//            hasFound = false;
//            while( !hasFound ) {
//                if( neighbor_grid_cell->type == 'w' ) {
//                    white_neighbor_cell = dynamic_cast<TransitionNode*>(neighbor_grid_cell);
//                    if( white_neighbor_cell->singleConnector == 0 ) {
//                        hasFound = true;
//                    }
//                    else {
//                        neighbor_grid_cell = white_neighbor_cell->quadrupleConnector1;
//                    }
//                }
//                else {
//                    break;
//                }
//            }
//
//        }
//    }
//        //Aqui, a célula vizinha tem um nivel de refinamento menor, entao eh mais simples.
//    else {
//        if(neighbor_grid_cell->level <= level && (neighbor_grid_cell->type == 'w') ) {
//            white_neighbor_cell = static_cast<TransitionNode*>(neighbor_grid_cell);
//            hasFound = false;
//            while( !hasFound ) {
//                if( neighbor_grid_cell->type == 'w' ) {
//                    white_neighbor_cell = dynamic_cast<TransitionNode*>(neighbor_grid_cell);
//                    if( white_neighbor_cell->singleConnector == 0 ) {
//                        hasFound = true;
//                    }
//                    else {
//                        neighbor_grid_cell = white_neighbor_cell->singleConnector;
//                    }
//                }
//                else {
//                    break;
//                }
//            }
//        }
//    }
//
//    //Tratamos somente os pontos interiores da malha.
//    if( ( neighbor_grid_cell->type == 'b' ) && ( neighbor_grid_cell->active == true ) )	{
//
//        black_neighbor_cell = static_cast<CellNode*>(neighbor_grid_cell);
//
//        if ( black_neighbor_cell->halfFaceLength < leastDistance )
//            leastDistance = black_neighbor_cell->halfFaceLength;
//
//        localFlux = ( v - black_neighbor_cell->v ) / ( 2 * leastDistance );
//
//        lock();
//
//        switch( direction ) {
//            case 's':
//                if ( localFlux > southFlux )
//                    southFlux += localFlux;
//                break;
//
//            case 'n':
//                if ( localFlux > northFlux )
//                    northFlux += localFlux;
//                break;
//
//            case 'e':
//                if ( localFlux > eastFlux )
//                    eastFlux += localFlux;
//                break;
//
//            case 'w':
//                if ( localFlux > westFlux )
//                    westFlux += localFlux;
//                break;
//
//            case 'f':
//                if ( localFlux > frontFlux )
//                    frontFlux += localFlux;
//                break;
//
//            case 'b':
//                if ( localFlux > backFlux )
//                    backFlux += localFlux;
//                break;
//        }
//
//        unlock();
//
//        black_neighbor_cell->lock();
//
//        switch( direction ) {
//            case 's':
//                if ( localFlux > black_neighbor_cell->northFlux )
//                    black_neighbor_cell->northFlux += localFlux;
//                break;
//
//            case 'n':
//                if ( localFlux > black_neighbor_cell->southFlux )
//                    black_neighbor_cell->southFlux += localFlux;
//                break;
//
//            case 'e':
//                if ( localFlux > black_neighbor_cell->westFlux )
//                    black_neighbor_cell->westFlux += localFlux;
//                break;
//
//            case 'w':
//                if ( localFlux > black_neighbor_cell->eastFlux )
//                    black_neighbor_cell->eastFlux += localFlux;
//                break;
//
//            case 'f':
//                if ( localFlux > black_neighbor_cell->backFlux )
//                    black_neighbor_cell->backFlux += localFlux;
//                break;
//
//            case 'b':
//                if ( localFlux > black_neighbor_cell->frontFlux )
//                    black_neighbor_cell->frontFlux += localFlux;
//                break;
//        }
//
//        black_neighbor_cell->unlock();
//
//
//    }
//}

