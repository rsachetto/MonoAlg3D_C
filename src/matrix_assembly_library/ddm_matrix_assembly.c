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
