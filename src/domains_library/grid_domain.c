//
// Created by sachetto on 01/10/17.
//

#include "domain_helpers.h"

#include "../utils/erros_helpers.h"
#include "../utils/logfile_utils.h"
#include "../main/config/domain_config.h"
#include "../libraries_common/config_helpers.h"
#include <assert.h>

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_mesh) {

    double start_h = config->start_h;

    int num_layers = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, num_layers, config->config_data.config, "num_layers");


    double side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, side_length, config->config_data.config, "side_length");


    double real_side_length = start_h * 2.0f;
    double max_h = start_h * num_layers;

    while (real_side_length < side_length) {
        real_side_length *= 2.0f;
    }

    print_to_stdout_and_file ("Initial cube side length: %lf µm x %lf µm x %lf µm\n", real_side_length,
                              real_side_length, real_side_length);
    print_to_stdout_and_file ("Loading plain mesh with %lf µm x %lf µm x %lf µm using dx %lf µm\n", side_length,
                              side_length, max_h, start_h);

    int num_steps = get_num_refinement_steps_to_discretization (real_side_length, start_h);

    initialize_and_construct_grid (the_grid, real_side_length, 7);

    if ((real_side_length / 2.0f) > max_h) {
        double aux = real_side_length / 2.0f;

        for (int i = 0; i < num_steps - 3; i++) {
            set_plain_domain (the_grid, real_side_length, real_side_length, aux);
            refine_grid (the_grid, 1);
            aux = aux / 2.0f;
        }

        refine_grid (the_grid, 3);

    } else {
        refine_grid (the_grid, num_steps);
    }

    set_plain_domain (the_grid, side_length, side_length, max_h);

    int i;
    for (i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells (the_grid);
    }
}

SET_SPATIAL_DOMAIN(initialize_grid_with_human_mesh) {

    config->start_h = 800.0;
    bool fibrotic = false;

    char *mesh_file;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(mesh_file, config->config_data.config, "mesh_file");

    char *fibrotic_char;
    GET_PARAMETER_VALUE_CHAR(fibrotic_char, config->config_data.config, "fibrotic");
    if (fibrotic_char != NULL) {
        fibrotic = ((strcmp (fibrotic_char, "yes") == 0) || (strcmp (fibrotic_char, "true") == 0));
    }

    double minx = -1;
    double maxx = -1;
    double miny = -1;
    double maxy = -1;
    double minz = -1;
    double maxz = -1;

    bool success;
    GET_PARAMETER_NUMERIC_VALUE(double, minx, config->config_data.config, "min_x", success);
    GET_PARAMETER_NUMERIC_VALUE(double, minx, config->config_data.config, "max_x", success);
    GET_PARAMETER_NUMERIC_VALUE(double, minx, config->config_data.config, "min_y", success);
    GET_PARAMETER_NUMERIC_VALUE(double, minx, config->config_data.config, "max_y", success);
    GET_PARAMETER_NUMERIC_VALUE(double, minx, config->config_data.config, "min_z", success);
    GET_PARAMETER_NUMERIC_VALUE(double, minx, config->config_data.config, "max_z", success);

    initialize_and_construct_grid (the_grid, 204800, 7);
    refine_grid (the_grid, 7);

    bool full_mesh = !((minx >= 0) && (maxx >= 0) && (miny >= 0) && (maxy >= 0) && (minz >= 0) && (maxz >= 0));
    if (full_mesh) {
        print_to_stdout_and_file ("Loading Human Heart Mesh\n");
        set_custom_mesh (the_grid, mesh_file, 2025252, fibrotic);
    } else {
        print_to_stdout_and_file ("Loading Human Heart Sub Mesh\n");
        set_custom_mesh_with_bounds (the_grid, mesh_file, 2025252, minx, maxx, miny, maxy, minz, maxz, fibrotic);
    }

    print_to_stdout_and_file ("Cleaning grid\n");
    int i;
    for (i = 0; i < 7; i++) {
        derefine_grid_inactive_cells (the_grid);
    }

    if (fibrotic) {
        refine_fibrotic_cells (the_grid);
        refine_fibrotic_cells (the_grid);
        refine_fibrotic_cells (the_grid);

        refine_border_zone_cells (the_grid);
        refine_border_zone_cells (the_grid);
        refine_border_zone_cells (the_grid);
        // set_human_mesh_fibrosis(the_grid, phi, scar_file, seed, scar_center_x, scar_center_y, scar_center_z);
    }

    free (mesh_file);
}

SET_SPATIAL_DOMAIN(initialize_grid_with_rabbit_mesh) {

    config->start_h = 250.0;

    char *mesh_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(mesh_file, config->config_data.config, "mesh_file");

    initialize_and_construct_grid (the_grid, 64000.0, 7);
    refine_grid (the_grid, 7);

    print_to_stdout_and_file ("Loading Rabbit Heart Mesh\n");

    set_custom_mesh (the_grid, mesh_file, 470197, false);

    print_to_stdout_and_file ("Cleaning grid\n");
    int i;
    for (i = 0; i < 6; i++) {
        derefine_grid_inactive_cells (the_grid);
    }
    free (mesh_file);
}

SET_SPATIAL_DOMAIN(initialize_grid_with_mouse_mesh) {

    char *mesh_file = NULL;

    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(mesh_file, config->config_data.config, "mesh_file");

    double start_h = config->start_h;

    assert (the_grid);

    initialize_and_construct_grid (the_grid, 6400.0, 7);

    refine_grid (the_grid, 5);

    print_to_stdout_and_file ("Loading Mouse Heart Mesh\n");

    set_custom_mesh (the_grid, mesh_file, 96195, false);

    print_to_stdout_and_file ("Cleaning grid\n");

    int i;
    for (i = 0; i < 5; i++) {
        derefine_grid_inactive_cells (the_grid);
    }

    if (start_h == 50.0) {
        print_to_stdout_and_file ("Refining Mesh to 50um\n");
        refine_grid (the_grid, 1);
    } else if (start_h == 25.0) {
        print_to_stdout_and_file ("Refining Mesh to 25um\n");
        refine_grid (the_grid, 1);
    }
}

SET_SPATIAL_DOMAIN(initialize_grid_with_benchmark_mesh) {

    double side_length;

    double start_h = config->start_h;

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

        while (grid_cell != 0) {
            if (grid_cell->active) {
                set_cell_not_changeable (grid_cell, start_h);
            }
            grid_cell = grid_cell->next;
        }
    }
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_fibrotic_mesh) {

    bool success;

    double phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, phi, config->config_data.config, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE(unsigned, seed, config->config_data.config, "seed", success);
    if(!success) seed = 0;

    initialize_grid_with_plain_mesh (the_grid, config);
    set_plain_fibrosis (the_grid, phi, seed);
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_and_sphere_fibrotic_mesh) {

    double phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, phi, config->config_data.config, "phi" );

    double plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, plain_center, config->config_data.config, "plain_center" );

    double sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, sphere_radius, config->config_data.config, "sphere_radius" );

    double border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, border_zone_size, config->config_data.config, "border_zone_size" );

    double border_zone_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, border_zone_radius, config->config_data.config, "border_zone_radius");

    bool success;
    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE(unsigned, seed, config->config_data.config, "seed", success);
    if (!success) {
        seed = 0;
    }

    initialize_grid_with_plain_mesh (the_grid, config);
    set_plain_sphere_fibrosis (the_grid, phi, plain_center, sphere_radius, border_zone_size, border_zone_radius, seed);
}
