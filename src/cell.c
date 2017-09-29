//
// Created by sachetto on 29/09/17.
//

#include "cell.h"


void init_basic_cell_data(struct basic_cell_data *data, char type,
                          uint8_t level, float center_x,float center_y,
                          float center_z) {

    data->type    = type;
    data->level   = level;
    data->center_x = center_x;
    data->center_y = center_y;
    data->center_z = center_z;

}

void init_basic_cell_data_with_default_values(struct basic_cell_data *data, char type) {

    data->type    = type;
    data->level   = 1;
    data->center_x = 0.0;
    data->center_y = 0.0;
    data->center_z = 0.0;

}

void init_cell_node(struct cell_node *cell_node, uint8_t init_ode) {

    init_basic_cell_data_with_default_values(&(cell_node->cell_data), 'b');
    cell_node->active = 1;

    cell_node->bunch_number = 0;
    cell_node->north       = 0;
    cell_node->south       = 0;
    cell_node->east        = 0;
    cell_node->west        = 0;
    cell_node->front       = 0;
    cell_node->back        = 0;
    cell_node->previous    = 0;
    cell_node->next        = 0;

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

    init_basic_cell_data_with_default_values(&(transition_node->cell_data), 'w');

    transition_node->single_connector     = 0;
    transition_node->quadruple_connector1 = 0;
    transition_node->quadruple_connector2 = 0;
    transition_node->quadruple_connector3 = 0;
    transition_node->quadruple_connector4 = 0;
    transition_node->direction            = 0;

}

