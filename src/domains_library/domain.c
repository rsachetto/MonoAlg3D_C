//
// Created by sachetto on 01/10/17.
//

#include "domain_helpers.h"

#include "../libraries_common/config_helpers.h"
#include "../monodomain/config/domain_config.h"
#include "../utils/logfile_utils.h"
#include <assert.h>
#include <time.h>

#ifdef _MSC_VER
    #include <process.h>
    #define getpid _getpid
#else
    #include <unistd.h>
#endif

SET_SPATIAL_DOMAIN (initialize_grid_with_plain_mesh) {

    double start_h = config->start_h;

    int num_layers = 0;
    bool s;
    GET_PARAMETER_NUMERIC_VALUE(int, num_layers, config->config_data.config, "num_layers", s);

    double side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, side_length, config->config_data.config, "side_length");

    double real_side_length = start_h * 2.0f;

    double max_h;

    if(!s || num_layers == 0) {
        max_h = side_length;
    }

    else {
        max_h = start_h * num_layers;
    }

    while (real_side_length < side_length) {
        real_side_length *= 2.0f;
    }

    print_to_stdout_and_file ("Initial cube side length: %lf µm x %lf µm x %lf µm\n", real_side_length,
                              real_side_length, real_side_length);
    print_to_stdout_and_file ("Loading plain mesh with %lf µm x %lf µm x %lf µm using dx %lf µm\n", side_length,
                              side_length, max_h, start_h);

    int num_steps = get_num_refinement_steps_to_discretization (real_side_length, start_h);

    initialize_and_construct_grid (the_grid, real_side_length);

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

SET_SPATIAL_DOMAIN (initialize_grid_with_human_mesh_with_two_scars) {

    config->start_h = 800.0;
    bool fibrotic = false;

    char *mesh_file;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (mesh_file, config->config_data.config, "mesh_file");

    char *fibrotic_char;
    GET_PARAMETER_VALUE_CHAR (fibrotic_char, config->config_data.config, "fibrotic");
    if (fibrotic_char != NULL) {
        fibrotic = ((strcmp (fibrotic_char, "yes") == 0) || (strcmp (fibrotic_char, "true") == 0));
    }

    initialize_and_construct_grid (the_grid, 204800);
    refine_grid (the_grid, 7);

    print_to_stdout_and_file ("Loading Human Heart Mesh\n");
    set_custom_mesh (the_grid, mesh_file, 2025252, fibrotic);

    print_to_stdout_and_file ("Cleaning grid\n");
    int i;
    for (i = 0; i < 7; i++) {
        derefine_grid_inactive_cells (the_grid);
    }

    if (fibrotic) {

        //Here we refine the scar cells
        refine_fibrotic_cells (the_grid); //400 um
        refine_fibrotic_cells (the_grid); //200 um
        refine_fibrotic_cells (the_grid); //100 um

        //and the border zone
        refine_border_zone_cells (the_grid);
        refine_border_zone_cells (the_grid);
        refine_border_zone_cells (the_grid);

        char *scar_file_big;
        GET_PARAMETER_VALUE_CHAR (scar_file_big, config->config_data.config, "big_scar_file");

        char *scar_file_small;
        GET_PARAMETER_VALUE_CHAR (scar_file_small, config->config_data.config, "small_scar_file");

        if(scar_file_big) {
            print_to_stdout_and_file("Loading fibrosis patterns from file %s\n", scar_file_big);
            set_human_mesh_fibrosis_from_file(the_grid, 'b', scar_file_big, 2172089);
        }

        if(scar_file_small) {
            print_to_stdout_and_file("Loading fibrosis patterns from file %s\n", scar_file_small);
            set_human_mesh_fibrosis_from_file(the_grid, 's', scar_file_small, 845051);
        }

        if(!(scar_file_big || scar_file_small)) {

            double small_scar_center_x = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, small_scar_center_x, config->config_data.config, "small_scar_center_x");

            double small_scar_center_y = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, small_scar_center_y, config->config_data.config, "small_scar_center_y");

            double small_scar_center_z = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, small_scar_center_z, config->config_data.config, "small_scar_center_z");

            double big_scar_center_x = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, big_scar_center_x, config->config_data.config, "big_scar_center_x");

            double big_scar_center_y = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, big_scar_center_y, config->config_data.config, "big_scar_center_y");

            double big_scar_center_z = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, big_scar_center_z, config->config_data.config, "big_scar_center_z");

            double phi = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, phi, config->config_data.config, "phi");

            unsigned seed = 0;
            bool seed_success ;
            GET_PARAMETER_NUMERIC_VALUE (unsigned, seed, config->config_data.config, "seed", seed_success);
            if (!seed_success)
                seed = 0;

            print_to_stdout_and_file("Setting random fibrosis pattern\n");
            set_human_mesh_fibrosis(the_grid, phi, seed, big_scar_center_x, big_scar_center_y, big_scar_center_z,
                                    small_scar_center_x, small_scar_center_y, small_scar_center_z);
        }

    }

    free (mesh_file);
}

SET_SPATIAL_DOMAIN(initialize_grid_with_scar_wedge) {
    char *mesh_file;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (mesh_file, config->config_data.config, "mesh_file");

    char *scar_size;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (scar_size, config->config_data.config, "scar_size");

    double phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, phi, config->config_data.config, "phi");

    unsigned fib_seed = 0;
    bool success;
    GET_PARAMETER_NUMERIC_VALUE(unsigned, fib_seed, config->config_data.config, "seed", success);

    if (!success)
        fib_seed = (unsigned)time (NULL) + getpid ();

    srand (fib_seed);

    config->start_h = 800.0;
    uint8_t size_code;

    initialize_and_construct_grid (the_grid, 204800);
    refine_grid (the_grid, 7);

    if(strcmp(scar_size, "big") == 0) {
        print_to_stdout_and_file("Loading Human Heart Edge with big scar\n");
        set_custom_mesh_with_bounds(the_grid, mesh_file, 2025252, 79100, 121000, 66700, 106000, 11200, 61400, true);
        size_code = 0;
    }
    else if(strcmp(scar_size, "small") == 0) {
        print_to_stdout_and_file("Loading Human Heart Edge with small scar\n");
        set_custom_mesh_with_bounds(the_grid, mesh_file, 2025252, 30400, 81600, 59200, 103000, 13600, 48000, true);
        size_code = 1;
    }
    else {
        printf("Function: initialize_grid_with_scar_edge, invalid scar size %s. Valid sizes are big or small. Exiting!\n", scar_size);
        exit(EXIT_FAILURE);
    }


    print_to_stdout_and_file ("Cleaning grid\n");
    int i;
    for (i = 0; i < 7; i++) {
        derefine_grid_inactive_cells (the_grid);
    }

    refine_fibrotic_cells (the_grid);
    refine_fibrotic_cells (the_grid);
    refine_fibrotic_cells (the_grid);

    refine_border_zone_cells (the_grid);
    refine_border_zone_cells (the_grid);
    refine_border_zone_cells (the_grid);

    double scar_center_x;
    double scar_center_y;
    double scar_center_z;

    ////Fibrosis configuration

    //BIG SCAR
    if(size_code == 0) {
        scar_center_x = 95300.0;
        scar_center_y = 81600.0;
        scar_center_z = 36800.0;
    }
    else {
        scar_center_x = 52469.0;
        scar_center_y = 83225.0;
        scar_center_z = 24791.0;
    }


    double bz_size = 0.0;
    double dist;

    print_to_stdout_and_file ("Using %u as seed\n", fib_seed);
    print_to_stdout_and_file("Calculating fibrosis using phi: %lf\n", phi);
    struct cell_node *grid_cell = the_grid->first_cell;
    
    while( grid_cell != 0 ) {

        if(grid_cell->active) {
            if(grid_cell->fibrotic) {
                grid_cell->can_change = false;
                double p = (double) (rand()) / (RAND_MAX);
                if (p < phi) grid_cell->active = false;
            }
            else if(grid_cell->border_zone) {
                double center_x = grid_cell->center_x;
                double center_y = grid_cell->center_y;
                double center_z = grid_cell->center_z;
                dist =  sqrt((center_x - scar_center_x)*(center_x - scar_center_x) + (center_y - scar_center_y)*(center_y - scar_center_y)  + (center_z - scar_center_z)*(center_z - scar_center_z)  );
                if(dist > bz_size) {
                    bz_size = dist;
                }
            }

        }
        grid_cell = grid_cell->next;
    }

    grid_cell = the_grid->first_cell;
    while( grid_cell != 0 ) {

        if(grid_cell->active) {
            if(grid_cell->border_zone) {
                double center_x = grid_cell->center_x;
                double center_y = grid_cell->center_y;
                double center_z = grid_cell->center_z;
                dist =  sqrt((center_x - scar_center_x)*(center_x - scar_center_x) + (center_y - scar_center_y)*(center_y - scar_center_y)  + (center_z - scar_center_z)*(center_z - scar_center_z)  );
                dist = dist/bz_size;

                double phi_local = phi - phi*dist;
                double p = (double) (rand()) / (RAND_MAX);

                if (p < phi_local) grid_cell->active = false;

                grid_cell->can_change = false;

            }

        }
        grid_cell = grid_cell->next;
    }

    free(mesh_file);
    free(scar_size);
}


SET_SPATIAL_DOMAIN (initialize_grid_with_rabbit_mesh) {

    config->start_h = 250.0;

    char *mesh_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (mesh_file, config->config_data.config, "mesh_file");

    initialize_and_construct_grid (the_grid, 64000.0);
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

SET_SPATIAL_DOMAIN (initialize_grid_with_mouse_mesh) {

    char *mesh_file = NULL;

    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (mesh_file, config->config_data.config, "mesh_file");

    double start_h = config->start_h;

    assert (the_grid);

    initialize_and_construct_grid (the_grid, 6400.0);

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

SET_SPATIAL_DOMAIN (initialize_grid_with_benchmark_mesh) {

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

    initialize_and_construct_grid (the_grid, side_length);
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

SET_SPATIAL_DOMAIN (initialize_grid_with_plain_fibrotic_mesh) {

    bool success;

    double phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, phi, config->config_data.config, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE (unsigned, seed, config->config_data.config, "seed", success);
    if (!success)
        seed = 0;

    initialize_grid_with_plain_mesh (config, the_grid);
    set_plain_fibrosis (the_grid, phi, seed);
}

SET_SPATIAL_DOMAIN (initialize_grid_with_plain_and_sphere_fibrotic_mesh) {

    double phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, phi, config->config_data.config, "phi");

    double plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, plain_center, config->config_data.config, "plain_center");

    double sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, sphere_radius, config->config_data.config, "sphere_radius");

    double border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, border_zone_size, config->config_data.config,
                                                 "border_zone_size");

    double border_zone_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (double, border_zone_radius, config->config_data.config,
                                                 "border_zone_radius");

    bool success;
    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE (unsigned, seed, config->config_data.config, "seed", success);
    if (!success) {
        seed = 0;
    }

    initialize_grid_with_plain_mesh (config, the_grid);
    set_plain_sphere_fibrosis (the_grid, phi, plain_center, sphere_radius, border_zone_size, border_zone_radius, seed);
}
