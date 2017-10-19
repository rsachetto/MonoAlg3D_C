//
// Created by sachetto on 01/10/17.
//

#include "grid_domain.h"
#include "domain_helpers.h"

#include "../utils/erros_helpers.h"
#include "../utils/logfile_utils.h"
#include <assert.h>
#include <string.h>


void initialize_grid_with_human_mesh (struct grid *the_grid, struct domain_config *domain_config) {

    domain_config->start_h = 800.0;
    bool fibrotic = false;

    char *mesh_file = string_hash_search(domain_config->config_data.config, "mesh_file");
    if(mesh_file == NULL) {
        report_parameter_error_on_function("initialize_grid_with_human_mesh", "mesh_file");
    }

    char *fibrotic_char = string_hash_search(domain_config->config_data.config, "fibrotic");
    if(fibrotic_char != NULL) {
        fibrotic = ((strcmp(fibrotic_char, "yes") == 0) || (strcmp(fibrotic_char, "true") == 0));
    }


    double minx = -1;
    double maxx = -1;
    double miny = -1;
    double maxy = -1;
    double minz = -1;
    double maxz = -1;

    char *config_char = string_hash_search(domain_config->config_data.config, "min_x");
    if(config_char)
        minx  = atof(config_char);
    free(config_char);

    config_char = string_hash_search(domain_config->config_data.config, "max_x");
    if(config_char)
        maxx  = atof(config_char);
    free(config_char);

    config_char = string_hash_search(domain_config->config_data.config, "min_y");
    if(config_char)
        miny  = atof(config_char);
    free(config_char);

    config_char = string_hash_search(domain_config->config_data.config, "max_y");
    if(config_char)
        maxy  = atof(config_char);
    free(config_char);

    config_char = string_hash_search(domain_config->config_data.config, "min_z");
    if(config_char)
        minz  = atof(config_char);
    free(config_char);

    config_char = string_hash_search(domain_config->config_data.config, "max_z");
    if(config_char)
        maxz  = atof(config_char);
    free(config_char);

    initialize_and_construct_grid (the_grid, 204800, 7);
    refine_grid (the_grid, 7);


    bool full_mesh = !((minx >= 0) && (maxx >= 0) && (miny >= 0) && (maxy >= 0) && (minz >= 0) && (maxz >= 0));
    if(full_mesh) {
        print_to_stdout_and_file ("Loading Human Heart Mesh\n");
        set_custom_mesh(the_grid, mesh_file, 2025252, fibrotic);
    }
    else {
        print_to_stdout_and_file ("Loading Human Heart Sub Mesh\n");
        set_custom_mesh_with_bounds (the_grid, mesh_file, 2025252, minx, maxx, miny, maxy, minz, maxz, fibrotic);
    }

    print_to_stdout_and_file ("Cleaning grid\n");
    int i;
    for (i = 0; i < 7; i++) {
        derefine_grid_inactive_cells (the_grid);
    }

    if(fibrotic) {
        refine_fibrotic_cells(the_grid);
        refine_fibrotic_cells(the_grid);
        refine_fibrotic_cells(the_grid);

        refine_border_zone_cells(the_grid);
        refine_border_zone_cells(the_grid);
        refine_border_zone_cells(the_grid);
        //set_human_mesh_fibrosis(the_grid, phi, scar_file, seed, scar_center_x, scar_center_y, scar_center_z);
    }
    free(mesh_file);



}

void initialize_grid_with_rabbit_mesh (struct grid *the_grid, struct domain_config *domain_config) {

    domain_config->start_h = 250.0;

    char *mesh_file = string_hash_search(domain_config->config_data.config, "mesh_file");
    if(mesh_file == NULL) {
        report_parameter_error_on_function("initialize_grid_with_rabbit_mesh", "mesh_file");
    }

    initialize_and_construct_grid (the_grid, 32000.0, 7);
    refine_grid (the_grid, 6);

    print_to_stdout_and_file ("Loading Rabbit Heart Mesh\n");

    set_custom_mesh(the_grid, mesh_file, 470197, false);

    print_to_stdout_and_file ("Cleaning grid\n");
    int i;
    for (i = 0; i < 6; i++) {
        derefine_grid_inactive_cells (the_grid);
    }
    free(mesh_file);
}

void initialize_grid_with_mouse_mesh (struct grid *the_grid, struct domain_config *domain_config) {

    char *mesh_file = string_hash_search(domain_config->config_data.config, "mesh_file");
    if(mesh_file == NULL) {
        report_parameter_error_on_function("initialize_grid_with_mouse_mesh", "mesh_file");
    }

    double start_h = domain_config->start_h;

    assert(the_grid);

    initialize_and_construct_grid (the_grid, 6400.0, 7);

    refine_grid (the_grid, 5);

    print_to_stdout_and_file ("Loading Mouse Heart Mesh\n");

    set_custom_mesh(the_grid, mesh_file, 96195, false);

    print_to_stdout_and_file ("Cleaning grid\n");

    int i;
    for (i = 0; i < 5; i++) {
        derefine_grid_inactive_cells (the_grid);
    }

    if(start_h == 50.0) {
        print_to_stdout_and_file("Refining Mesh to 50um\n");
        refine_grid(the_grid, 1);
    }
    else if (start_h == 25.0) {
        print_to_stdout_and_file("Refining Mesh to 25um\n");
        refine_grid(the_grid, 1);
    }

}



void initialize_grid_with_benchmark_mesh (struct grid *the_grid, struct domain_config *domain_config) {

    double side_length;

    double start_h = domain_config->start_h;

    print_to_stdout_and_file ("Loading N-Version benchmark mesh using dx %lf um\n", start_h);
    if ((start_h == 100.0) || (start_h == 200.0)) {
        side_length = 25600.0;
    } else if (start_h == 125.0 || start_h == 250.0 || start_h == 500.0) {
        side_length = 32000.0;
    } else {
        fprintf (stderr, "initialize_grid_with_benchmark_mesh: invalid value of start_h (initial "
                "discretization). Exiting!");
        exit (10);
    }

    initialize_and_construct_grid (the_grid, side_length, 7);
    int num_steps = get_num_refinement_steps_to_discretization (side_length, start_h);

    refine_grid (the_grid, num_steps);
    set_benchmark_domain (the_grid);

    print_to_stdout_and_file ("Cleaning grid\n");
    int i;

    for (i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells (the_grid);
    }

    if (the_grid->adaptive) {
        struct cell_node *grid_cell;
        grid_cell = the_grid->first_cell;

        while( grid_cell != 0 ) {
            if(grid_cell->active) {
                set_cell_not_changeable(grid_cell, start_h);
            }
            grid_cell = grid_cell->next;
        }
    }


}

void initialize_grid_with_plain_mesh (struct grid *the_grid, struct domain_config* config) {

    double start_h = config->start_h;

    char *config_char = string_hash_search(config->config_data.config, "num_layers");
    if(config_char == NULL) {
        report_parameter_error_on_function("initialize_grid_with_benchmark_mesh", "num_layers");
    }
    int num_layers = atoi(config_char);
    free(config_char);

    config_char = string_hash_search(config->config_data.config, "side_length");
    if(config_char == NULL) {
        report_parameter_error_on_function("initialize_grid_with_benchmark_mesh", "side_length");
    }
    double desired_side_length = atof(config_char);
    free(config_char);


    double real_side_length = start_h * 2.0f;
    double max_h = start_h * num_layers;

    while (real_side_length < desired_side_length) {
        real_side_length *= 2.0f;
    }

    print_to_stdout_and_file("Initial cube side length: %lf µm x %lf µm x %lf µm\n", real_side_length, real_side_length, real_side_length);
    print_to_stdout_and_file("Loading plain mesh with %lf µm x %lf µm x %lf µm using dx %lf µm\n", desired_side_length, desired_side_length, max_h, start_h);

    int num_steps = get_num_refinement_steps_to_discretization(real_side_length, start_h);

    initialize_and_construct_grid(the_grid, real_side_length, 7);

    if ((real_side_length / 2.0f) > max_h) {
        double aux = real_side_length/2.0f;

        for (int i = 0; i < num_steps - 3; i++) {
            set_plain_domain(the_grid, real_side_length, real_side_length, aux);
            refine_grid(the_grid, 1);
            aux = aux / 2.0f;
        }

        refine_grid(the_grid, 3);


    } else {
        refine_grid(the_grid, num_steps);
    }

    set_plain_domain(the_grid, desired_side_length, desired_side_length, max_h);

    int i;
    for (i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells(the_grid);
    }


}

void initialize_grid_with_plain_fibrotic_mesh(struct grid *the_grid, struct domain_config* config) {

    char *config_char = string_hash_search(config->config_data.config, "phi");
    if(config_char == NULL) {
        report_parameter_error_on_function("initialize_grid_with_plain_fibrotic_mesh", "phi");
    }
    double phi = atof(config_char);
    free(config_char);

    config_char = string_hash_search(config->config_data.config, "seed");
    unsigned seed = 0;
    if(config_char) {
        seed = (unsigned )atoi(config_char);
        free(config_char);
    }


    initialize_grid_with_plain_mesh(the_grid, config);
    set_plain_fibrosis(the_grid, phi, seed);

}


void initialize_grid_with_plain_and_sphere_fibrotic_mesh(struct grid *the_grid, struct domain_config* config) {

    char *config_char = string_hash_search(config->config_data.config, "phi");
    if(config_char == NULL) {
        report_parameter_error_on_function("initialize_grid_with_plain_and_sphere_fibrotic_mesh", "phi");
    }

    double phi = atof(config_char);
    free(config_char);

    config_char = string_hash_search(config->config_data.config, "plain_center");
    if(config_char == NULL) {
        report_parameter_error_on_function("initialize_grid_with_plain_and_sphere_fibrotic_mesh", "plain_center");
    }
    double plain_center = atof(config_char);
    free(config_char);

    config_char = string_hash_search(config->config_data.config, "sphere_radius");
    if(config_char == NULL) {
        report_parameter_error_on_function("initialize_grid_with_plain_and_sphere_fibrotic_mesh", "sphere_radius");
    }
    double sphere_radius = atof(config_char);
    free(config_char);

    config_char = string_hash_search(config->config_data.config, "border_zone_size");
    if(config_char == NULL) {
        report_parameter_error_on_function("initialize_grid_with_plain_and_sphere_fibrotic_mesh", "border_zone_size");
    }
    double border_zone_size = atof(config_char);
    free(config_char);


    config_char = string_hash_search(config->config_data.config, "border_zone_radius");
    if(config_char == NULL) {
        report_parameter_error_on_function("initialize_grid_with_plain_and_sphere_fibrotic_mesh", "border_zone_radius");
    }
    double border_zone_radius = atof(config_char);
    free(config_char);

    config_char = string_hash_search(config->config_data.config, "seed");
    unsigned seed = 0;
    if(config_char) {
        seed = (unsigned )atoi(config_char);
        free(config_char);
    }


    initialize_grid_with_plain_mesh(the_grid, config);
    set_plain_sphere_fibrosis(the_grid, phi, plain_center,sphere_radius,border_zone_size, border_zone_radius, seed);


}

