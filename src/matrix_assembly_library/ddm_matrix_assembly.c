//
// Created by bergolho on 25/07/19.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../libraries_common/config_helpers.h"
#include "../monodomain/constants.h"
#include "../utils/utils.h"
#include "../single_file_libraries/stb_ds.h"

static struct element fill_element_ddm (uint32_t position, char direction, real_cpu dx, real_cpu dy, real_cpu dz,\
                                 const real_cpu sigma_x, const real_cpu sigma_y, const real_cpu sigma_z,\
                                 const real_cpu kappa_x, const real_cpu kappa_y, const real_cpu kappa_z,\
                                 const real_cpu dt,\
                                struct element *cell_elements);

void create_sigma_low_block(struct cell_node* ac,real_cpu x_left, real_cpu x_right, real_cpu y_down, real_cpu y_up,double b_sigma_x, double b_sigma_y, double b_sigma_z,double sigma_factor)
{
    real_cpu x = ac->center_x;
    real_cpu y = ac->center_y;

    
    if( (x<=x_right) && (x>= x_left) && (y>= y_down) && (y<= y_up))
    {
            ac->sigma_x = b_sigma_x*sigma_factor;
            ac->sigma_y = b_sigma_y*sigma_factor;
            ac->sigma_z = b_sigma_z*sigma_factor;
    }
	
}

void calculate_kappa_elements(struct monodomain_solver *the_solver, struct grid *the_grid,\
                    const real_cpu cell_length_x, const real_cpu cell_length_y, const real_cpu cell_length_z)
{
	uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    
    #pragma omp parallel for
    for (uint32_t i = 0; i < num_active_cells; i++) 
	{
        ac[i]->kappa_x = KAPPA(beta,cm,cell_length_x,ac[i]->dx);
        ac[i]->kappa_y = KAPPA(beta,cm,cell_length_y,ac[i]->dy);
        ac[i]->kappa_z = KAPPA(beta,cm,cell_length_z,ac[i]->dz);
	}
	
}

struct element fill_element_ddm (uint32_t position, char direction, real_cpu dx, real_cpu dy, real_cpu dz,\
                                 const real_cpu sigma_x, const real_cpu sigma_y, const real_cpu sigma_z,\
                                 const real_cpu kappa_x, const real_cpu kappa_y, const real_cpu kappa_z,\
                                 const real_cpu dt,\
                                struct element *cell_elements) 
{

    real_cpu multiplier;

    struct element new_element;
    new_element.column = position;
    new_element.direction = direction;

    // Z direction
    if(direction == 'n') 
    { 
        multiplier = ((dx * dy) / dz);
        new_element.value = ( multiplier * (-sigma_z - (kappa_z / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_z + (kappa_z / dt)) );
    } 
    // Z direction
    else if(direction == 's') 
    { 
        multiplier = ((dx * dy) / dz);
        new_element.value = ( multiplier * (-sigma_z - (kappa_z / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_z + (kappa_z / dt)) );
    } 
    // Y direction
    else if(direction == 'e') 
    { 
        multiplier = ((dx * dz) / dy);
        new_element.value = ( multiplier * (-sigma_y - (kappa_y / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_y + (kappa_y / dt)) );
    } 
    // Y direction
    else if(direction == 'w') 
    { 
        multiplier = ((dx * dz) / dy);
        new_element.value = ( multiplier * (-sigma_y - (kappa_y / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_y + (kappa_y / dt)) );
    }
    // X direction 
    else if(direction == 'f') 
    { 
        multiplier = ((dy * dz) / dx);
        new_element.value = ( multiplier * (-sigma_x - (kappa_x / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_x + (kappa_x / dt)) );
    } 
    // X direction
    else if(direction == 'b') 
    { 
        multiplier = ((dy * dz) / dx);
        new_element.value = ( multiplier * (-sigma_x - (kappa_x / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_x + (kappa_x / dt)) );
    }
    
    return new_element;
}

static void fill_discretization_matrix_elements_ddm (struct cell_node *grid_cell, void *neighbour_grid_cell, real_cpu dt ,char direction)                                                 

{
	uint32_t position;
    bool has_found;
    real_cpu dx, dy, dz;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if(neighbour_grid_cell_level > grid_cell->cell_data.level) 
    {
        if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) 
        {
            has_found = false;
            while(!has_found) 
            {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) 
                {
                    white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                    if(white_neighbor_cell->single_connector == NULL) 
                    {
                        has_found = true;
                    } 
                    else 
                    {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } 
                else 
                {
                    break;
                }
            }
        }
    } 
    else 
    {
        if(neighbour_grid_cell_level <= grid_cell->cell_data.level &&
           (neighbour_grid_cell_type == TRANSITION_NODE_TYPE)) 
           {
            has_found = false;
            while(!has_found) 
            {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) 
                {
                    white_neighbor_cell = (struct transition_node *)(neighbour_grid_cell);
                    if(white_neighbor_cell->single_connector == 0) 
                    {
                        has_found = true;
                    } 
                    else 
                    {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } 
                else 
                {
                    break;
                }
            }
        }
    }

    // We care only with the interior points
    if(neighbour_grid_cell_type == CELL_NODE_TYPE) 
    {

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if(black_neighbor_cell->active) 
        {

            real_cpu sigma_x1 = grid_cell->sigma_x;
            real_cpu sigma_x2 = black_neighbor_cell->sigma_x;
            real_cpu sigma_x = 0.0;
            
            if(sigma_x1 != 0.0 && sigma_x2 != 0.0) 
            {
                sigma_x = (2.0f * sigma_x1 * sigma_x2) / (sigma_x1 + sigma_x2);
            }

            real_cpu sigma_y1 = grid_cell->sigma_y;
            real_cpu sigma_y2 = black_neighbor_cell->sigma_y;
            real_cpu sigma_y = 0.0;

            if(sigma_y1 != 0.0 && sigma_y2 != 0.0) 
            {
                sigma_y = (2.0f * sigma_y1 * sigma_y2) / (sigma_y1 + sigma_y2);
            }

            real_cpu sigma_z1 = grid_cell->sigma_z;
            real_cpu sigma_z2 = black_neighbor_cell->sigma_z;
            real_cpu sigma_z = 0.0;

            if(sigma_z1 != 0.0 && sigma_z2 != 0.0) 
            {
                sigma_z = (2.0f * sigma_z1 * sigma_z2) / (sigma_z1 + sigma_z2);
            }
            
            if(black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) 
            {
                dx = black_neighbor_cell->dx;
                dy = black_neighbor_cell->dy;
                dz = black_neighbor_cell->dz;
            }
            else 
            {
                dx = grid_cell->dx;
                dy = grid_cell->dy;
                dz = grid_cell->dz;
            }

            lock_cell_node(grid_cell);

            struct element *cell_elements = grid_cell->elements;
            position = black_neighbor_cell->grid_position;

            size_t max_elements = arrlen(cell_elements);
            bool insert = true;

            for(size_t i = 1; i < max_elements; i++) 
            {
                if(cell_elements[i].column == position) 
                {
                    insert = false;
                    break;
                }
            }

            if(insert) 
            {

                struct element new_element = fill_element_ddm(position, direction,\
																dx, dy, dz,\
                                                                sigma_x, sigma_y, sigma_z,\
																grid_cell->kappa_x,grid_cell->kappa_y,grid_cell->kappa_z,\
																dt, cell_elements);

                new_element.cell = black_neighbor_cell;
                arrput(grid_cell->elements, new_element);
            }
            unlock_cell_node(grid_cell);

            lock_cell_node(black_neighbor_cell);
            cell_elements = black_neighbor_cell->elements;
            position = grid_cell->grid_position;

            max_elements = arrlen(cell_elements);

            insert = true;
            for(size_t i = 1; i < max_elements; i++) 
            {
                if(cell_elements[i].column == position) 
                {
                    insert = false;
                    break;
                }
            }

            if(insert) 
            {

                struct element new_element = fill_element_ddm(position, direction,\
																dx, dy, dz,\
                                                                sigma_x, sigma_y, sigma_z,\
																grid_cell->kappa_x,grid_cell->kappa_y,grid_cell->kappa_z,\
																dt, cell_elements);

                new_element.cell = grid_cell;
                arrput(black_neighbor_cell->elements, new_element);
            }

            unlock_cell_node(black_neighbor_cell);
        }
    }
}

void initialize_diagonal_elements(struct monodomain_solver *the_solver, struct grid *the_grid) 
{

    real_cpu alpha, dx, dy, dz;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    int i;

    #pragma omp parallel for private(alpha, dx, dy, dz)
    for(i = 0; i < num_active_cells; i++) 
    {
        dx = ac[i]->dx;
        dy = ac[i]->dy;
        dz = ac[i]->dz;

        alpha = ALPHA(beta, cm, dt, dx, dy, dz);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if(ac[i]->elements)
            arrfree(ac[i]->elements);

        ac[i]->elements = NULL;

        arrsetcap(ac[i]->elements, 7);
        arrput(ac[i]->elements, element);
    }
}

int randRange(int n) {
    int limit;
    int r;

    limit = RAND_MAX - (RAND_MAX % n);

    while((r = rand()) >= limit)
        ;

    return r % n;
}

// This function will read the fibrotic regions and for each cell that is inside the region and we will
// reduce its conductivity value based on the 'sigma_factor'.
ASSEMBLY_MATRIX (heterogenous_fibrotic_sigma_with_factor_ddm_assembly_matrix)
{
    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    char *fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file, config->config_data.config, "fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config->config_data.config, "size");	

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");

    real cell_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_x, config->config_data.config, "cell_length_x");

    real cell_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_y, config->config_data.config, "cell_length_y");

    real cell_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_z, config->config_data.config, "cell_length_z");
      
    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data.config, "sigma_factor");    

    // Calculate the kappa values on each cell of th grid
	calculate_kappa_elements(the_solver,the_grid,cell_length_x,cell_length_y,cell_length_z);
  
    if(!sigma_initialized) 
    {
        #pragma omp parallel for
        for (uint32_t i = 0; i < num_active_cells; i++) 
        {
            ac[i]->sigma_x = sigma_x;
            ac[i]->sigma_y = sigma_y;
            ac[i]->sigma_z = sigma_z;
        }

        sigma_initialized = true;
    }

    // Read and store the fibrosis locations
    FILE *file = fopen(fib_file, "r");

    if(!file) 
    {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) 
    {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 7);
        if(scar_mesh[i] == NULL) 
        {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    uint32_t i = 0;
    while (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3], &scar_mesh[i][4], &scar_mesh[i][5], &scar_mesh[i][6]) != EOF)
    {
        i++;
    }

    fclose(file); 

    uint32_t num_fibrotic_regions = i;

    // Pass through all the cells of the grid and check if its center is inside the current
    // fibrotic region
    #pragma omp parallel for
    for(int j = 0; j < num_fibrotic_regions; j++) 
    {

        struct cell_node *grid_cell = the_grid->first_cell;
    
        real_cpu b_center_x = scar_mesh[j][0];
        real_cpu b_center_y = scar_mesh[j][1];

        real_cpu b_h_dx = scar_mesh[j][3];
        real_cpu b_h_dy = scar_mesh[j][4];

        bool active = (bool) (scar_mesh[j][6]);
       
        while(grid_cell != 0) 
        {
            if (grid_cell->active)
            {
                real_cpu center_x = grid_cell->center_x;
                real_cpu center_y = grid_cell->center_y;
                real_cpu half_dy = grid_cell->dy/2.0;
                real_cpu half_dx = grid_cell->dx/2.0;

                struct point_3d p;
                struct point_3d q;
                
                p.x = b_center_y + b_h_dy;
                p.y = b_center_y - b_h_dy; 
                
                q.x = b_center_x + b_h_dx;
                q.y = b_center_x - b_h_dx; 

                // Check if the current cell is inside the fibrotic region
                if (center_x > q.y && center_x < q.x && center_y > p.y && center_y < p.x)  
                {
                    if(active == 0)
                    {
                        grid_cell->sigma_x = sigma_x * sigma_factor;
                        grid_cell->sigma_y = sigma_y * sigma_factor;
                        grid_cell->sigma_z = sigma_z * sigma_factor;
                    }
                }
            }    
            grid_cell = grid_cell->next;
        }    
    }

    printf("[!] Using DDM formulation\n");
    printf("[X] Cell length = %.10lf || sigma_x = %.10lf || dx = %.10lf || kappa_x = %.10lf\n",\
            cell_length_x,ac[0]->sigma_x,ac[0]->dx,ac[0]->kappa_x);
    printf("[Y] Cell length = %.10lf || sigma_y = %.10lf || dy = %.10lf || kappa_y = %.10lf\n",\
            cell_length_y,ac[0]->sigma_y,ac[0]->dy,ac[0]->kappa_y);
    printf("[Z] Cell length = %.10lf || sigma_z = %.10lf || dz = %.10lf || kappa_z = %.10lf\n",\
            cell_length_z,ac[0]->sigma_z,ac[0]->dz,ac[0]->kappa_z); 


    #pragma omp parallel for
    for(int i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->south, the_solver->dt,'s');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->north, the_solver->dt,'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->east, the_solver->dt,'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->west, the_solver->dt,'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->front, the_solver->dt,'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->back, the_solver->dt,'b');
    }


    for(int k = 0; k < fib_size; k++) 
    {
        free(scar_mesh[k]);
    }

    free(scar_mesh);
}

ASSEMBLY_MATRIX(homogenous_ddm_assembly_matrix) 
{
    static bool sigma_initialized = false;

    struct cell_node *grid_cell;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");

    real cell_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_x, config->config_data.config, "cell_length_x");

    real cell_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_y, config->config_data.config, "cell_length_y");

    real cell_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_z, config->config_data.config, "cell_length_z");

    // Calculate the kappa values on each cell of th grid
	calculate_kappa_elements(the_solver,the_grid,cell_length_x,cell_length_y,cell_length_z);

    printf("[!] Using DDM formulation\n");
    printf("[X] Cell length = %.10lf || sigma_x = %.10lf || dx = %.10lf || kappa_x = %.10lf\n",\
            cell_length_x,ac[0]->sigma_x,ac[0]->dx,ac[0]->kappa_x);
    printf("[Y] Cell length = %.10lf || sigma_y = %.10lf || dy = %.10lf || kappa_y = %.10lf\n",\
            cell_length_y,ac[0]->sigma_y,ac[0]->dy,ac[0]->kappa_y);
    printf("[Z] Cell length = %.10lf || sigma_z = %.10lf || dz = %.10lf || kappa_z = %.10lf\n",\
            cell_length_z,ac[0]->sigma_z,ac[0]->dz,ac[0]->kappa_z);

    // Initialize the conductivities of each cell
    if (!sigma_initialized)
    {
	    grid_cell = the_grid->first_cell;
	    while(grid_cell != 0) 
	    {

		    if(grid_cell->active) 
		    {
	    		grid_cell->sigma_x = sigma_x;
	    		grid_cell->sigma_y = sigma_y;
	    		grid_cell->sigma_z = sigma_z;
		    }
		    grid_cell = grid_cell->next;
    	}

	    sigma_initialized = true;
    }

    #pragma omp parallel for
    for(int i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->south, the_solver->dt,'s');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->north, the_solver->dt,'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->east, the_solver->dt,'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->west, the_solver->dt,'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->front, the_solver->dt,'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->back, the_solver->dt,'b');
    }
    
}

ASSEMBLY_MATRIX(sigma_low_region_triangle_ddm_tiny) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");
    
    real cell_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_x, config->config_data.config, "cell_length_x");

    real cell_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_y, config->config_data.config, "cell_length_y");

    real cell_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_z, config->config_data.config, "cell_length_z");    
    

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data.config, "sigma_factor");

    real sigma_factor_2 = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor_2, config->config_data.config, "sigma_factor_2");

    real side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, side_length, config->config_data.config, "side_length");

    
    // Calculate the kappas for the DDM
    calculate_kappa_elements(the_solver,the_grid,cell_length_x,cell_length_y,cell_length_z);

    //~ bool inside;

	// Initialize the conductivities
    #pragma omp parallel for
	for (i = 0; i < num_active_cells; i++) 
    {
		ac[i]->sigma_x = sigma_x;
		ac[i]->sigma_y = sigma_y;
		ac[i]->sigma_z = sigma_z;
	}


    printf("[!] Using DDM formulation\n");
    printf("[X] Cell length = %.10lf || sigma_x = %.10lf || dx = %.10lf || kappa_x = %.10lf\n",\
            cell_length_x,ac[0]->sigma_x,ac[0]->dx,ac[0]->kappa_x);
    printf("[Y] Cell length = %.10lf || sigma_y = %.10lf || dy = %.10lf || kappa_y = %.10lf\n",\
            cell_length_y,ac[0]->sigma_y,ac[0]->dy,ac[0]->kappa_y);
    printf("[Z] Cell length = %.10lf || sigma_z = %.10lf || dz = %.10lf || kappa_z = %.10lf\n",\
            cell_length_z,ac[0]->sigma_z,ac[0]->dz,ac[0]->kappa_z);

//regiao ao redor do meio com sigma menor//
	
/*Aqui comeca a brincadeira de criar os triangulos
 * 
 * Pensei em dividir em 3 partes a primeira uma barra
 * 
 *    1 2 3
 *    *******
 *	  ******
 *    *****
 *    ****
 *    **
 * 
 *
 * 1 fazendo o canal
 * 
 *  	**
 *  	**
 *  	**
 *  	**
 * 		**
 * 
 * 2 considerar a abertura de maneira que pegue 5 celulas
 * 
 * 		*
 * 		*
 * 		*
 * 		*
 * 		
 * 3 pegar o triangulo
 * 
 * 		****
 * 		***
 * 		**
 * 		*
 * 
 * Dai (1)+(2)+(3) eh a area que desejamos
 * 
 * */


// set the square scar mesh...


	//~ real X_left  = side_length/3.0;
	//~ real X_right = (side_length/3.0+side_length/3.0);
	//~ real Y_down  = side_length/3.0;
	//~ real Y_up 	= (side_length/3.0+side_length/3.0);	
	
	//esquerda
	
	//~ real X_left  = side_length/6.0;
	//~ real X_right = (side_length/3.0+side_length/3.0);
	real X_left  = side_length/12.0;
	real X_right = (side_length/6.0+side_length/6.0);
	//~ real Y_down  = side_length/6.0;
	real Y_down  = 0.0;
	//~ real Y_up 	= (side_length/6.0+side_length/6.0+side_length/12.0);	
	real Y_up 	= (side_length/6.0+side_length/6.0+side_length/6.0);	
	

//General square slow sigma *aquela divisao por 3 pode mudar*...
	
	#pragma omp parallel for
	for (i = 0; i < num_active_cells; i++) 
	{
			double x_prev = ac[i]->center_x;
			double y_prev = ac[i]->center_y;
		
			if ((x_prev>= X_left ) && (x_prev <= X_right) && (y_prev>= Y_down) && (y_prev<= Y_up))
			{
				ac[i]->sigma_x = sigma_x * sigma_factor;
				ac[i]->sigma_y = sigma_y * sigma_factor;
				ac[i]->sigma_z = sigma_z * sigma_factor;
			}
	}	
	
//Vou recuparando os triangulos do square	
	#pragma omp parallel for
	for (i = 0; i < num_active_cells; i++) 
	{	

		double x = ac[i]->center_x;
		double y = ac[i]->center_y;	
		
		//Part 1
		int mid_scar_Y = (Y_down + Y_up)/12.0;
		int Part1_x_right = X_left + 10.*(cell_length_x);
		int Part2_x_left = Part1_x_right + cell_length_x;
		
		
		if( (y<=(mid_scar_Y+(cell_length_y/2.))) && (y> (mid_scar_Y- (cell_length_y/2.))) && (x>= X_left) && (x<= X_right))
		
		{
			ac[i]-> sigma_x = sigma_x;
			ac[i]-> sigma_y = sigma_y;
			ac[i]-> sigma_z = sigma_z;
		}
		
		
		//Part 1.5
		if( (y<=(mid_scar_Y+(cell_length_y))) && (y> (mid_scar_Y- 2.0*(cell_length_y))) && (x>= X_left) && (x<= X_left+5.0*cell_length_x))  //isto mudei
		{
			ac[i]-> sigma_x = sigma_x;
			ac[i]-> sigma_y = sigma_y;
			ac[i]-> sigma_z = sigma_z;
		}
		
		
		//Part 2
		int Part2_y_up = mid_scar_Y + (3.*cell_length_y);	//original
		int Part2_y_down = mid_scar_Y - (2.*cell_length_y); //original
		//~ int Part2_y_up = mid_scar_Y + (20.*cell_length_y);//modifiquei
		//~ int Part2_y_down = mid_scar_Y - (20.*cell_length_y);//modifiquei //talvez colocar o de cima
		
				
		if( (x>= Part2_x_left) && (x<= X_right) && (y<= Part2_y_up ) &&(y> Part2_y_down))
		{
			ac[i]-> sigma_x = sigma_x;
			ac[i]-> sigma_y = sigma_y;
			ac[i]-> sigma_z = sigma_z;
		}
		
		
		//Part 3 up, m coeficiente angular
		
		int Part3_x_left = Part2_x_left + cell_length_x;
		int Part3_y_up = Part2_y_up;
		real m_up = (Y_up-Part3_y_up)/(X_right-Part3_x_left);
		int Part3_x_left_2 = (Part3_x_left + X_right)/2.0;

		
		if((x>=Part3_x_left) && (x <= X_right) && (y <= (m_up*(x-Part3_x_left)+Part3_y_up)) && (y>= Part3_y_up) )
			{	
				
				if(x>=Part3_x_left_2)
				{
					ac[i]-> sigma_x = sigma_x*sigma_factor_2;
					ac[i]-> sigma_y = sigma_y*sigma_factor_2;
					ac[i]-> sigma_z = sigma_z*sigma_factor_2;
				}
				
				else
				{
					ac[i]-> sigma_x = sigma_x;
					ac[i]-> sigma_y = sigma_y;
					ac[i]-> sigma_z = sigma_z;	
				}
				
			}
		
		//Part 3 down, m coeficiente angular
		
		int Part3_y_down = Part2_y_down;
		
		real m_down = (Y_down - Part3_y_down)/(X_right-Part3_x_left);
		
		if((x>=Part3_x_left) && (x <= X_right) && (y >= (m_down*(x-Part3_x_left)+ Part3_y_down)) && (y<= Part3_y_down) && (y>=Y_down)) 
			{
				if(x>=Part3_x_left_2)
				{
					ac[i]-> sigma_x = sigma_x*sigma_factor_2;
					ac[i]-> sigma_y = sigma_y*sigma_factor_2;
					ac[i]-> sigma_z = sigma_z*sigma_factor_2;
				}
				
				else
				{
					ac[i]-> sigma_x = sigma_x;
					ac[i]-> sigma_y = sigma_y;
					ac[i]-> sigma_z = sigma_z;	
				}				
				
			}	
	
        //middle canal
		if ((y<=Part3_y_up) && (y>=Part3_y_down) && (x>=Part3_x_left_2) && (x<=X_right))
        {	
                ac[i]-> sigma_x = sigma_x*sigma_factor_2;
                ac[i]-> sigma_y = sigma_y*sigma_factor_2;
                ac[i]-> sigma_z = sigma_z*sigma_factor_2;
        }
			
	}

// Aqui termina a construcao dos triangulos	
	
    	

    #pragma omp parallel for
    for (i = 0; i < num_active_cells; i++)
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->south, the_solver->dt,'s');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->north, the_solver->dt,'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->east, the_solver->dt,'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->west, the_solver->dt,'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->front, the_solver->dt,'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->back, the_solver->dt,'b');
    }
}


ASSEMBLY_MATRIX(write_sigma_low_region_triangle_ddm_tiny) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
	struct cell_node *grid_cell;


    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");
    
    real cell_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_x, config->config_data.config, "cell_length_x");

    real cell_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_y, config->config_data.config, "cell_length_y");

    real cell_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_z, config->config_data.config, "cell_length_z");    
    

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data.config, "sigma_factor");

    real sigma_factor_2 = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor_2, config->config_data.config, "sigma_factor_2");

    real side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, side_length, config->config_data.config, "side_length");

    char *new_fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(new_fib_file, config->config_data.config, "rescaled_fibrosis_file");    

    //~ real scar_length = 0.0;
    //~ GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, scar_length, config->config_data.config, "scar_length");

	//~ calculate_kappa_elements(the_grid,cell_length_x,cell_length_y,cell_length_z);
	//~ calculate_kappa_elements(the_solver,the_grid,cell_length_x,cell_length_y,cell_length_z);



    //~ bool inside;

	#pragma omp parallel for
	for (i = 0; i < num_active_cells; i++) {
		ac[i]->sigma_x = sigma_x;
		ac[i]->sigma_y = sigma_y;
		ac[i]->sigma_z = sigma_z;
	}


    printf("[!] Using DDM formulation\n");
    printf("[X] Cell length = %.10lf || sigma_x = %.10lf || dx = %.10lf || kappa_x = %.10lf\n",\
            cell_length_x,ac[0]->sigma_x,ac[0]->dx,ac[0]->kappa_x);
    printf("[Y] Cell length = %.10lf || sigma_y = %.10lf || dy = %.10lf || kappa_y = %.10lf\n",\
            cell_length_y,ac[0]->sigma_y,ac[0]->dy,ac[0]->kappa_y);
    printf("[Z] Cell length = %.10lf || sigma_z = %.10lf || dz = %.10lf || kappa_z = %.10lf\n",\
            cell_length_z,ac[0]->sigma_z,ac[0]->dz,ac[0]->kappa_z);
//regiao ao redor do meio com sigma menor//
	
/*Aqui comeca a brincadeira de criar os triangulos
 * 
 * Pensei em dividir em 3 partes a primeira uma barra
 * 
 *    1 2 3
 *    *******
 *	  ******
 *    *****
 *    ****
 *    **
 * 
 * 
 * 1 fazendo o canal
 * 
 *  	**
 *  	**
 *  	**
 *  	**
 * 		**
 * 
 * 2 considerar a abertura de maneira que pegue 5 celulas
 * 
 * 		*
 * 		*
 * 		*
 * 		*
 * 		
 * 3 pegar o triangulo
 * 
 * 		****
 * 		***
 * 		**
 * 		*
 * 
 * Dai (1)+(2)+(3) eh a ahrea que desejamos
 * 
 * */


// set the square scar mesh...


	//~ real X_left  = side_length/3.0;
	//~ real X_right = (side_length/3.0+side_length/3.0);
	//~ real Y_down  = side_length/3.0;
	//~ real Y_up 	= (side_length/3.0+side_length/3.0);	
	
	//esquerda
	
	//~ real X_left  = side_length/6.0;
	//~ real X_right = (side_length/3.0+side_length/3.0);
	real X_left  = side_length/12.0;
	real X_right = (side_length/6.0+side_length/6.0);
	//~ real Y_down  = side_length/6.0;
	real Y_down  = 0.0;
	//~ real Y_up 	= (side_length/6.0+side_length/6.0+side_length/12.0);	
	real Y_up 	= (side_length/6.0+side_length/6.0+side_length/6.0);	


//General square slow sigma *aquela divisao por 3 pode mudar*...
	
	#pragma omp parallel for
	for (i = 0; i < num_active_cells; i++) 
	{
			double x_prev = ac[i]->center_x;
			double y_prev = ac[i]->center_y;
		
			if ((x_prev>= X_left ) && (x_prev <= X_right) && (y_prev>= Y_down) && (y_prev<= Y_up))
			{
				ac[i]->sigma_x = sigma_x * sigma_factor;
				ac[i]->sigma_y = sigma_y * sigma_factor;
				ac[i]->sigma_z = sigma_z * sigma_factor;
			}
	}	
	
//Vou recuparando os triangulos do square	
	
	#pragma omp parallel for
	for (i = 0; i < num_active_cells; i++) 
	{	

		double x = ac[i]->center_x;
		double y = ac[i]->center_y;	
		
		//Part 1
		

		int mid_scar_Y = (Y_down + Y_up)/12.0;
		int Part1_x_right = X_left + 10.*(cell_length_x);
		int Part2_x_left = Part1_x_right + cell_length_x;
		
		
		if( (y<=(mid_scar_Y+(cell_length_y/2.))) && (y> (mid_scar_Y- (cell_length_y/2.))) && (x>= X_left) && (x<= X_right))
		
		{
			ac[i]-> sigma_x = sigma_x;
			ac[i]-> sigma_y = sigma_y;
			ac[i]-> sigma_z = sigma_z;
		}
		
		
		//Part 1.5
		
		
		if( (y<=(mid_scar_Y+(cell_length_y))) && (y> (mid_scar_Y- 2.0*(cell_length_y))) && (x>= X_left) && (x<= X_left+5.0*cell_length_x))  //isto mudei
		{
			ac[i]-> sigma_x = sigma_x;
			ac[i]-> sigma_y = sigma_y;
			ac[i]-> sigma_z = sigma_z;
		}
		
		
		//~ //Part 2
		
		int Part2_y_up = mid_scar_Y + (3.*cell_length_y);	//original
		int Part2_y_down = mid_scar_Y - (2.*cell_length_y); //original
		//~ int Part2_y_up = mid_scar_Y + (20.*cell_length_y);//modifiquei
		//~ int Part2_y_down = mid_scar_Y - (20.*cell_length_y);//modifiquei //talvez colocar o de cima
		
				
		if( (x>= Part2_x_left) && (x<= X_right) && (y<= Part2_y_up ) &&(y> Part2_y_down))
		{
			ac[i]-> sigma_x = sigma_x;
			ac[i]-> sigma_y = sigma_y;
			ac[i]-> sigma_z = sigma_z;
		}
		
		
		//~ //Part 3 up, m coeficiente angular
		
		int Part3_x_left = Part2_x_left + cell_length_x;
		int Part3_y_up = Part2_y_up;
		real m_up = (Y_up-Part3_y_up)/(X_right-Part3_x_left);
		int Part3_x_left_2 = (Part3_x_left + X_right)/2.0;

		
		if((x>=Part3_x_left) && (x <= X_right) && (y <= (m_up*(x-Part3_x_left)+Part3_y_up)) && (y>= Part3_y_up) )
			{	
				
				if(x>=Part3_x_left_2)
				{
					ac[i]-> sigma_x = sigma_x*sigma_factor_2;
					ac[i]-> sigma_y = sigma_y*sigma_factor_2;
					ac[i]-> sigma_z = sigma_z*sigma_factor_2;
				}
				
				else
				{
					ac[i]-> sigma_x = sigma_x;
					ac[i]-> sigma_y = sigma_y;
					ac[i]-> sigma_z = sigma_z;	
				}
				
			}
		
		//~ //Part 3 down, m coeficiente angular
		
		int Part3_y_down = Part2_y_down;
		
		real m_down = (Y_down - Part3_y_down)/(X_right-Part3_x_left);
		
		if((x>=Part3_x_left) && (x <= X_right) && (y >= (m_down*(x-Part3_x_left)+ Part3_y_down)) && (y<= Part3_y_down) && (y>=Y_down)) 
			{
				if(x>=Part3_x_left_2)
				{
					ac[i]-> sigma_x = sigma_x*sigma_factor_2;
					ac[i]-> sigma_y = sigma_y*sigma_factor_2;
					ac[i]-> sigma_z = sigma_z*sigma_factor_2;
				}
				
				else
				{
					ac[i]-> sigma_x = sigma_x;
					ac[i]-> sigma_y = sigma_y;
					ac[i]-> sigma_z = sigma_z;	
				}				
				
			}	
	
	
	
	//middle canal

	
	
		if ((y<=Part3_y_up) && (y>=Part3_y_down) && (x>=Part3_x_left_2) && (x<=X_right))
			{	
					ac[i]-> sigma_x = sigma_x*sigma_factor_2;
					ac[i]-> sigma_y = sigma_y*sigma_factor_2;
					ac[i]-> sigma_z = sigma_z*sigma_factor_2;
			}
			
	}
		
    // Write the new grid configuration on the rescaled_fibrosis file
	FILE *fileW = fopen(new_fib_file, "w+");

    grid_cell = the_grid->first_cell;

    while(grid_cell != 0) 
    {

        if(grid_cell->active) 
        {
            // We reescale the cell position using the 'rescale_factor'
            double center_x = grid_cell->center_x ;
            double center_y = grid_cell->center_y ;
            double center_z = grid_cell->center_z ;
            double dx = grid_cell->dx;
            double dy = grid_cell->dy;
            double dz = grid_cell->dz;
            double w_sigma_x = grid_cell->sigma_x;
            double w_sigma_y = grid_cell->sigma_y;
            double w_sigma_z = grid_cell->sigma_z;
                
            // Then, we write only the fibrotic regions to the output file
            fprintf(fileW,"%g,%g,%g,%g,%g,%g,%g,%g,%g\n",center_x,center_y,center_z,dx/2.0,dy/2.0,dz/2.0,w_sigma_x,w_sigma_y,w_sigma_z);
            
        }
        grid_cell = grid_cell->next;
    }

	fclose(fileW);  	
		
    // We just leave the program after this ...
    print_to_stdout_and_file("[!] Finish writing fibrotic region file '%s'!\n",new_fib_file);
    exit(EXIT_SUCCESS);

}

ASSEMBLY_MATRIX (heterogenous_fibrotic_sigma_with_factor_ddm_assembly_matrix_add_sigma_reverse)
{
    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    char *fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file, config->config_data.config, "fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config->config_data.config, "size");	

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");

    real cell_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_x, config->config_data.config, "cell_length_x");

    real cell_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_y, config->config_data.config, "cell_length_y");

    real cell_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_z, config->config_data.config, "cell_length_z");
      
    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data.config, "sigma_factor");    

    // Calculate the kappa values on each cell of th grid
	calculate_kappa_elements(the_solver,the_grid,cell_length_x,cell_length_y,cell_length_z);
  
    if(!sigma_initialized) 
    {
        #pragma omp parallel for
        for (uint32_t i = 0; i < num_active_cells; i++) 
        {
            ac[i]->sigma_x = sigma_x;
            ac[i]->sigma_y = sigma_y;
            ac[i]->sigma_z = sigma_z;
        }

        sigma_initialized = true;
    }

    // Read and store the fibrosis locations
    FILE *file = fopen(fib_file, "r");

    if(!file) 
    {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) 
    {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 9);
        if(scar_mesh[i] == NULL) 
        {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    uint32_t i = 0;
    while (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3],\
																 &scar_mesh[i][4], &scar_mesh[i][5], &scar_mesh[i][6], &scar_mesh[i][7],\
																 &scar_mesh[i][8]) != EOF)
    {
        i++;
    }

    fclose(file); 

    uint32_t num_fibrotic_regions = i;

    // Pass through all the cells of the grid and check if its center is inside the current
    // fibrotic region
 
    struct cell_node *grid_cell = the_grid->first_cell;

    real aux_sigma_load_x;
    real aux_sigma_load_y;
    real aux_sigma_load_z;

    while(grid_cell != 0)
    {
		if (grid_cell->active)
		{				
			real_cpu center_x = grid_cell-> center_x;
			real_cpu center_y = grid_cell-> center_y;
				
			real_cpu half_dx = grid_cell->dx/2.0;
			real_cpu half_dy = grid_cell->dy/2.0;
			aux_sigma_load_x = grid_cell->sigma_x;
			aux_sigma_load_y = grid_cell->sigma_y;
			aux_sigma_load_z = grid_cell->sigma_z;

            struct point_3d p;
            struct point_3d q;

            p.x = center_y + half_dy;
            p.y = center_y - half_dy;

            q.x = center_x + half_dx;
            q.y = center_x - half_dx;            
			
			for (int j = 0; j < num_fibrotic_regions; j++)
				{
        
					real_cpu b_center_x = scar_mesh[j][0];
					real_cpu b_center_y = scar_mesh[j][1];
					real_cpu b_sigma_x  = scar_mesh[j][6];
					real_cpu b_sigma_y  = scar_mesh[j][7];
					real_cpu b_sigma_z  = scar_mesh[j][8];

							
					if (b_center_x > q.y && b_center_x < q.x && b_center_y > p.y && b_center_y < p.x)
					{
                        
	
                            aux_sigma_load_x = b_sigma_x;
                            aux_sigma_load_y = b_sigma_y;
                            aux_sigma_load_z = b_sigma_z;
                            
                            //armazenar os valores dos sigmas ao redor e utilizar a meia harmonica para resolvers
	
					}
				
				
				
				}
				
            //~ grid_cell->sigma_x = 2.0*((sigma_x * aux_sigma_load_x)/(sigma_x + aux_sigma_load_x));
            //~ grid_cell->sigma_y = 2.0*((sigma_y * aux_sigma_load_y)/(sigma_y + aux_sigma_load_y));
            //~ grid_cell->sigma_z = 2.0*((sigma_z * aux_sigma_load_z)/(sigma_z + aux_sigma_load_z));
		   grid_cell->sigma_x = aux_sigma_load_x;
            grid_cell->sigma_y = aux_sigma_load_y;
            grid_cell->sigma_z = aux_sigma_load_z;


		}


        grid_cell = grid_cell->next;

	}

    printf("[!] Using DDM formulation\n");
    printf("[X] Cell length = %.10lf || sigma_x = %.10lf || dx = %.10lf || kappa_x = %.10lf\n",\
            cell_length_x,ac[0]->sigma_x,ac[0]->dx,ac[0]->kappa_x);
    printf("[Y] Cell length = %.10lf || sigma_y = %.10lf || dy = %.10lf || kappa_y = %.10lf\n",\
            cell_length_y,ac[0]->sigma_y,ac[0]->dy,ac[0]->kappa_y);
    printf("[Z] Cell length = %.10lf || sigma_z = %.10lf || dz = %.10lf || kappa_z = %.10lf\n",\
            cell_length_z,ac[0]->sigma_z,ac[0]->dz,ac[0]->kappa_z); 


    #pragma omp parallel for
    for(int i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->south, the_solver->dt,'s');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->north, the_solver->dt,'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->east, the_solver->dt,'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->west, the_solver->dt,'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->front, the_solver->dt,'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements_ddm(ac[i], ac[i]->back, the_solver->dt,'b');
    }


    for(int k = 0; k < fib_size; k++) 
    {
        free(scar_mesh[k]);
    }

    free(scar_mesh);
}

ASSEMBLY_MATRIX(heterogenous_fibrotic_region_file_write_using_seed)
{

    static bool sigma_initialized = false;
    int num;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    struct cell_node *grid_cell;

    initialize_diagonal_elements(the_solver, the_grid);

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data.config, "sigma_factor");    
    
    real sigma_factor_2 = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor_2, config->config_data.config, "sigma_factor_2");        

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data.config, "phi");
    
    real_cpu phi_2 = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi_2, config->config_data.config, "phi_2");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config->config_data.config, "seed");

    char *new_fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(new_fib_file, config->config_data.config, "rescaled_fibrosis_file");

    double rescale_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real,rescale_factor, config->config_data.config, "rescale_factor");
    
    real x_shift = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, x_shift, config->config_data.config, "x_shift");        
    
    real y_shift = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, y_shift, config->config_data.config, "y_shift");            

    // Write the new fibrotic region file
	FILE *fileW = fopen(new_fib_file, "w+");

    grid_cell = the_grid->first_cell;

    // Initialize the random the generator with the same seed used by the original model
    srand(seed);
    while(grid_cell != 0) 
    {

        if(grid_cell->active) 
        {
            real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
            if(p < phi) 
            {
                // We reescale the cell position using the 'rescale_factor'
				    grid_cell->sigma_x = sigma_x * sigma_factor_2;
				    grid_cell->sigma_y = sigma_y * sigma_factor_2;
				    grid_cell->sigma_z = sigma_z * sigma_factor_2;

                
                
                // Then, we write only the fibrotic regions to the output file
             
            }
            
            else
            {
				if((p>=phi)&&(p<=phi_2))
				{
				    grid_cell->sigma_x = sigma_x * sigma_factor;
				    grid_cell->sigma_y = sigma_y * sigma_factor;
				    grid_cell->sigma_z = sigma_z * sigma_factor;
				    
				}
				
				else
				{
				    grid_cell->sigma_x = sigma_x;
				    grid_cell->sigma_y = sigma_y;
				    grid_cell->sigma_z = sigma_z;		
				}	
				
			}
				
			double center_x = grid_cell->center_x;
			double center_y = grid_cell->center_y; 
			double center_z = grid_cell->center_z; 

			double dx = grid_cell->dx ;
			double dy = grid_cell->dy ;
			double dz = grid_cell->dz ;
			
			center_x = center_x + x_shift;
			center_y = center_y + y_shift;
			
			fprintf(fileW,"%g,%g,%g,%g,%g,%g,%g,%g,%g\n",center_x,center_y,center_z,dx/2.0,dy/2.0,dz/2.0,grid_cell->sigma_x,grid_cell->sigma_y,grid_cell->sigma_z);

        }
        grid_cell = grid_cell->next;
    }

    fclose(fileW);  	
		
    // We just leave the program after this ...
    print_to_stdout_and_file("[!] Finish writing fibrotic region file '%s'!\n",new_fib_file);
    exit(EXIT_SUCCESS);
}

ASSEMBLY_MATRIX(sigma_low_region_triangle_ddm_tiny_random_write) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    char *fib_file_1 = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file_1, config->config_data.config, "fibrosis_file_1");

    char *fib_file_2 = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file_2, config->config_data.config, "fibrosis_file_2");

    char *fib_file_3 = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file_3, config->config_data.config, "fibrosis_file_3");
    
    char *new_fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(new_fib_file, config->config_data.config, "new_fib_file");    

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config->config_data.config, "size");	    

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");
    
    real cell_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_x, config->config_data.config, "cell_length_x");

    real cell_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_y, config->config_data.config, "cell_length_y");

    real cell_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_z, config->config_data.config, "cell_length_z");    
    
    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data.config, "sigma_factor");

    real sigma_factor_2 = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor_2, config->config_data.config, "sigma_factor_2");

    real side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, side_length, config->config_data.config, "side_length");

	calculate_kappa_elements(the_solver,the_grid,cell_length_x,cell_length_y,cell_length_z);

	#pragma omp parallel for
	for (i = 0; i < num_active_cells; i++) 
	{
		ac[i]->sigma_x = sigma_x;
		ac[i]->sigma_y = sigma_y;
		ac[i]->sigma_z = sigma_z;
	}


    printf("[!] Using DDM formulation\n");
    printf("[X] Cell length = %.10lf || sigma_x = %.10lf || dx = %.10lf || kappa_x = %.10lf\n",\
            cell_length_x,ac[0]->sigma_x,ac[0]->dx,ac[0]->kappa_x);
    printf("[Y] Cell length = %.10lf || sigma_y = %.10lf || dy = %.10lf || kappa_y = %.10lf\n",\
            cell_length_y,ac[0]->sigma_y,ac[0]->dy,ac[0]->kappa_y);
    printf("[Z] Cell length = %.10lf || sigma_z = %.10lf || dz = %.10lf || kappa_z = %.10lf\n",\
            cell_length_z,ac[0]->sigma_z,ac[0]->dz,ac[0]->kappa_z);

    // Allocate memory to store the grid configuration
    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) 
    {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 9);
        if(scar_mesh[i] == NULL) 
        {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    FILE *file;
    i = 0;

    // Read the FIRST fibrosis file
    file = fopen(fib_file_1, "r");

    if(!file) 
    {
        printf("Error opening file %s!!\n", fib_file_1);
        exit(0);
    }

    
    while (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3],\
																 &scar_mesh[i][4], &scar_mesh[i][5], &scar_mesh[i][6], &scar_mesh[i][7],\
																 &scar_mesh[i][8]) != EOF)
    {
        i++;
    }

    fclose(file); 

    // Read the SECOND fibrosis file
    file = fopen(fib_file_2, "r");

    if(!file) 
    {
        printf("Error opening file %s!!\n", fib_file_2);
        exit(0);
    }

    
    while (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3],\
																 &scar_mesh[i][4], &scar_mesh[i][5], &scar_mesh[i][6], &scar_mesh[i][7],\
																 &scar_mesh[i][8]) != EOF)
    {
        i++;
    }

    fclose(file); 

    // Read the THIRD fibrosis file
    file = fopen(fib_file_3, "r");

    if(!file) 
    {
        printf("Error opening file %s!!\n", fib_file_3);
        exit(0);
    }

    
    while (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3],\
																 &scar_mesh[i][4], &scar_mesh[i][5], &scar_mesh[i][6], &scar_mesh[i][7],\
																 &scar_mesh[i][8]) != EOF)
    {
        i++;
    }

    fclose(file); 

    uint32_t num_fibrotic_regions = i;

    // Update the cells that are inside of the scar regions
   #pragma omp parallel for
    for(int j = 0; j < num_fibrotic_regions; j++) 
    {

        struct cell_node *grid_cell = the_grid->first_cell;
    
        real_cpu b_center_x = scar_mesh[j][0];
        real_cpu b_center_y = scar_mesh[j][1];

        real_cpu b_h_dx = scar_mesh[j][3];
        real_cpu b_h_dy = scar_mesh[j][4];
		real aux_sigma_load_x = scar_mesh[j][6];
		real aux_sigma_load_y  = scar_mesh[j][7];
		       
        while(grid_cell != 0) 
        {
            if (grid_cell->active)
            {
                real_cpu center_x = grid_cell->center_x;
                real_cpu center_y = grid_cell->center_y;
                real_cpu half_dy = grid_cell->dy/2.0;
                real_cpu half_dx = grid_cell->dx/2.0;

                struct point_3d p;
                struct point_3d q;
                
                p.x = b_center_y + b_h_dy;
                p.y = b_center_y - b_h_dy; 
                
                q.x = b_center_x + b_h_dx;
                q.y = b_center_x - b_h_dx; 

                // Check if the current cell is inside the fibrotic region
                if (center_x > q.y && center_x < q.x && center_y > p.y && center_y < p.x)  
                {
                        grid_cell->sigma_x = aux_sigma_load_x;
                        grid_cell->sigma_y = aux_sigma_load_y;
                }
            }    
            grid_cell = grid_cell->next;
        }    
    }
 
    // Until here, we have a grid with both the column and horizontal regions ...

    // Start building the channel block
	real X_left  = side_length/12.0;
	real X_left_1  = X_left + side_length/12.0;
	real X_right = 3000;
	real Y_down  = 0.0;
	real Y_up 	= 1200;
	
	#pragma omp parallel for
	for (i = 0; i < num_active_cells; i++) 
	{	

		double x = ac[i]->center_x;
		double y = ac[i]->center_y;	
		
    // Region 1 starts here ...

		//Part 1 - (Channel construction)

		int mid_scar_Y = (Y_down + 7500)/12.0;
		int X_right_1 = X_left + cell_length_x;
		int Part1_x_right = X_left + 10.*(cell_length_x);
		int Part2_x_left = Part1_x_right + cell_length_x;
		
		// Decrease the conductivity of every cell inside the region
		create_sigma_low_block(ac[i],X_left,Part1_x_right,Y_down,Y_up,sigma_x,sigma_y,sigma_z,sigma_factor);
	
		// Create the little channel
		create_sigma_low_block(ac[i],X_left,Part1_x_right,mid_scar_Y - cell_length_y/2.,mid_scar_Y+(cell_length_y/4.), sigma_x,sigma_y,sigma_z,1.0 );
		create_sigma_low_block(ac[i],X_left,X_left+2.0*cell_length_x,mid_scar_Y- 3.0*(cell_length_y),mid_scar_Y, sigma_x,sigma_y,sigma_z,1.0);
	
		//Part 2 - (Channel opening)
		
		int Part2_y_up = mid_scar_Y + (3.*cell_length_y);	//original
		int Part2_y_down = mid_scar_Y - (4.*cell_length_y); //original
		
		int X_left_2 = X_left+5.0*cell_length_x;
		int X_right_2 = X_left+10.0*cell_length_x;
		
		
		create_sigma_low_block(ac[i], X_left_2,X_right_2,Part2_y_down,Part2_y_up, sigma_x,sigma_y,sigma_z,1.0);		
				
	
		//Part 3 - triangle opening

        // First step
		int X_left_3 = X_left_2 + cell_length_x;
		int X_right_3 = X_right_2;
		int Y_up_3 = Part2_y_up + 5*cell_length_y;
		int Y_down_3 = Part2_y_down - 5*cell_length_y;
		real_cpu m_up = (Y_up_3-Part2_y_up)/(X_right_3-X_left_3);


		create_sigma_low_block(ac[i], X_left_3,X_right_3,Y_down_3,Y_up_3,sigma_x,sigma_y,sigma_z,1.0);		
			
        
		// Second step
		int X_left_4 = X_left_3 + cell_length_x;
		int X_right_4 = X_right_3;
		int Y_up_4 = Y_up_3 + 5*cell_length_y;
		int Y_down_4 = Y_down_3 - 5*cell_length_y;		
		
		create_sigma_low_block(ac[i], X_left_4,X_right_4,Y_down_4,Y_up_4,sigma_x,sigma_y,sigma_z,1.0);		

		
        // Third step
		int X_left_5 = X_left_4 + cell_length_x;
		int X_right_5 = X_right_4;
		int Y_up_5 = Y_up_4 + 5*cell_length_y;
		int Y_down_5 = Y_down_4 - 5*cell_length_y;		
		
		create_sigma_low_block(ac[i], X_left_5,X_right_5,Y_down_5,Y_up_5,sigma_x,sigma_y,sigma_z,1.0);	


		// Fourth step
		int X_left_6 = X_left_5 + cell_length_x;
		int X_right_6 = X_right_5;
		int Y_up_6 = Y_up_5 + 5*cell_length_y;
		int Y_down_6 = Y_down_5 - 5*cell_length_y;		
		
		create_sigma_low_block(ac[i], X_left_6,X_right_6,Y_down_6,Y_up_6,sigma_x,sigma_y,sigma_z,1.0);	
		

	// Region 2 starts here ...
		
        // The horizontal barrier
		int X_left_7 = 1900;
		int X_right_7 = 3900;
		int Y_down_7 = 1500;
		int Y_up_7 = 1900;
		
		create_sigma_low_block(ac[i], X_left_7,X_right_7,Y_down_7,Y_up_7,sigma_x,sigma_y,sigma_z,sigma_factor_2);	

        // The vertical barrier
		int X_left_8 = 3500;
		int X_right_8 = 3900;
		int Y_down_8 = 0;
		int Y_up_8 = 1700;
		
		create_sigma_low_block(ac[i], X_left_8,X_right_8,Y_down_8,Y_up_8,sigma_x,sigma_y,sigma_z,sigma_factor_2);	
		

	// Region 3 starts here ...
		
        // First upper block
		int X_left_9 = 1900;
		int X_right_9 = 2300;
		int Y_down_9 = 7900;
		int Y_up_9 = 8300;
		
		create_sigma_low_block(ac[i], X_left_9,X_right_9,Y_down_9,Y_up_9,sigma_x,sigma_y,sigma_z,sigma_factor_2);	

        // Second upper block
		int X_left_9_1 = 2600;
		int X_right_9_1 = 3000;
		int Y_down_9_1 = 7900;
		int Y_up_9_1 = 8300;
		
		create_sigma_low_block(ac[i], X_left_9_1,X_right_9_1,Y_down_9_1,Y_up_9_1,sigma_x,sigma_y,sigma_z,sigma_factor_2);	
		
        // Third upper block (reset)
		int X_left_10 = 2600;
		int X_right_10 = 3000;
		int Y_down_10 = 8600;
		int Y_up_10 = 9000;
		
		create_sigma_low_block(ac[i], X_left_10,X_right_10,Y_down_10,Y_up_10,sigma_x,sigma_y,sigma_z,1.0);			
		
        // Fourth upper block
		int X_left_11 = 3800;
		int X_right_11 = 4200;
		int Y_down_11 = 8200;
		int Y_up_11 = 8600;
		
		create_sigma_low_block(ac[i], X_left_11,X_right_11,Y_down_11,Y_up_11,sigma_x,sigma_y,sigma_z,sigma_factor_2);					
		
	}

    // Write the grid configuration to an output file
	FILE *fileW = fopen(new_fib_file, "w+");

    struct cell_node *grid_cell = the_grid->first_cell;

    // Initialize the random the generator with the same seed used by the original model
    while(grid_cell != 0) 
    {

        if(grid_cell->active) 
        {
                double center_x = grid_cell->center_x ;
                double center_y = grid_cell->center_y ;
                double center_z = grid_cell->center_z ;
                double dx = grid_cell->dx;
                double dy = grid_cell->dy;
                double dz = grid_cell->dz;
                double w_sigma_x = grid_cell->sigma_x;
                double w_sigma_y = grid_cell->sigma_y;
                double w_sigma_z = grid_cell->sigma_z;
                
                fprintf(fileW,"%g,%g,%g,%g,%g,%g,%g,%g,%g\n",center_x,center_y,center_z,dx/2.0,dy/2.0,dz/2.0,w_sigma_x,w_sigma_y,w_sigma_z);
        }
        grid_cell = grid_cell->next;
    }
	
	
	
	fclose(fileW);  	
		
    // We just leave the program after this ...
    print_to_stdout_and_file("[!] Finish writing fibrotic region file '%s'!\n",new_fib_file);
    exit(EXIT_SUCCESS);

}
