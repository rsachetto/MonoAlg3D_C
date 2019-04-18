//
// Created by sachetto on 01/10/17.
//

#include "domain_helpers.h"

#include "../config/domain_config.h"
#include "../libraries_common/common_data_structures.h"
#include "../libraries_common/config_helpers.h"
#include "../string/sds.h"
#include "../utils/file_utils.h"
#include "../single_file_libraries/stb_ds.h"
#include <assert.h>
#include <time.h>

#ifdef _MSC_VER
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

SET_SPATIAL_DOMAIN(initialize_grid_with_cuboid_mesh) {

    real_cpu start_dx = config->start_dx;
    real_cpu start_dy = config->start_dy;
    real_cpu start_dz = config->start_dz;

    real_cpu side_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length_x, config->config_data.config, "side_length_x");

    real_cpu side_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length_y, config->config_data.config, "side_length_y");

    real_cpu side_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length_z, config->config_data.config, "side_length_z");

    real_cpu real_side_length_x;
    real_cpu real_side_length_y;
    real_cpu real_side_length_z;

    int success =
        calculate_cuboid_side_lengths(start_dx, start_dy, start_dz, side_length_x, side_length_y, side_length_z,
                                      &real_side_length_x, &real_side_length_y, &real_side_length_z);

    if(!success) {
        return 0;
    }

    print_to_stdout_and_file("Initial mesh side length: %lf µm x %lf µm x %lf µm\n", real_side_length_x,
                             real_side_length_y, real_side_length_z);
    print_to_stdout_and_file(
        "Loading cuboid mesh with %lf µm x %lf µm x %lf µm using dx %lf µm, dy %lf µm, dz %lf µm\n", side_length_x,
        side_length_y, side_length_z, start_dx, start_dy, start_dz);

    int num_steps = get_num_refinement_steps_to_discretization(real_side_length_x, start_dx);

    initialize_and_construct_grid(the_grid, real_side_length_x, real_side_length_y, real_side_length_z);

    //TODO: maybe to this in all directions????
    if((real_side_length_z / 2.0f) > side_length_z) {
        real_cpu aux = real_side_length_z / 2.0f;
        int remaining_refinements = num_steps;

        for(int i = 0; i < num_steps; i++) {
            set_cuboid_domain(the_grid, real_side_length_x, real_side_length_y, aux);
            refine_grid(the_grid, 1);
            derefine_grid_inactive_cells(the_grid);
            aux = aux / 2.0f;
            remaining_refinements--;
            if(aux <= side_length_z) break;
        }

        refine_grid(the_grid, remaining_refinements);

    } else {
        refine_grid(the_grid, num_steps);
    }

    set_cuboid_domain(the_grid, side_length_x, side_length_y, side_length_z);

    int i;
    for(i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    translate_visible_mesh_to_origin(the_grid);

    return 1;
}

SET_SPATIAL_DOMAIN (initialize_grid_with_square_mesh) {

    int num_layers = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, num_layers, config->config_data.config, "num_layers");

    real_cpu side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real_cpu, side_length, config->config_data.config, "side_length");

    sds sx_char = sdscatprintf(sdsempty(), "%lf", side_length);
    sds sy_char = sdscatprintf(sdsempty(), "%lf", side_length);
    sds sz_char = sdscatprintf(sdsempty(), "%lf", config->start_dz*num_layers);


    shput(config->config_data.config, "side_length_x", strdup(sx_char));
    shput(config->config_data.config, "side_length_y", strdup(sy_char));
    shput(config->config_data.config, "side_length_z", strdup(sz_char));

    return initialize_grid_with_cuboid_mesh(config, the_grid);


}

SET_SPATIAL_DOMAIN(initialize_grid_with_cable_mesh) {

    real_cpu start_dx = config->start_dx;
    real_cpu start_dy = config->start_dy;
    real_cpu start_dz = config->start_dz;

    real_cpu cable_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, cable_length, config->config_data.config, "cable_length");

    real_cpu real_side_length_x;
    real_cpu real_side_length_y;
    real_cpu real_side_length_z;

    int success = calculate_cuboid_side_lengths(start_dx, start_dy, start_dz, cable_length, start_dy, start_dz,
                                                &real_side_length_x, &real_side_length_y, &real_side_length_z);

    if(!success) {
        return 0;
    }

    print_to_stdout_and_file("Loading cable mesh with %lf µm using dx %lf µm, dy %lf µm, dz %lf µm\n", cable_length,
                             start_dy, start_dz);

    int num_steps = get_num_refinement_steps_to_discretization(real_side_length_x, start_dx);

    initialize_and_construct_grid(the_grid, real_side_length_x, real_side_length_y, real_side_length_z);

    refine_grid(the_grid, num_steps);

    set_cuboid_domain(the_grid, cable_length, start_dy, start_dz);

    int i;
    for(i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    translate_visible_mesh_to_origin(the_grid);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_human_mesh_with_two_scars) {

    config->start_dx = 800.0;
    config->start_dy = 800.0;
    config->start_dx = 800.0;

    bool fibrotic = false;

    char *mesh_file;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(mesh_file, config->config_data.config, "mesh_file");

    GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(fibrotic, config->config_data.config, "fibrotic");

    initialize_and_construct_grid(the_grid, 204800, 204800, 204800);
    refine_grid(the_grid, 7);

    print_to_stdout_and_file("Loading Human Heart Mesh\n");
    set_custom_mesh(the_grid, mesh_file, 2025252, fibrotic);

    print_to_stdout_and_file("Cleaning grid\n");
    int i;
    for(i = 0; i < 7; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    if(fibrotic) {

        // Here we refine the scar cells
        refine_fibrotic_cells(the_grid); // 400 um
        refine_fibrotic_cells(the_grid); // 200 um
        refine_fibrotic_cells(the_grid); // 100 um

        // and the border zone
        refine_border_zone_cells(the_grid);
        refine_border_zone_cells(the_grid);
        refine_border_zone_cells(the_grid);

        char *scar_file_big = NULL;
        GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(scar_file_big, config->config_data.config, "big_scar_file");

        char *scar_file_small = NULL;
        GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(scar_file_small, config->config_data.config, "small_scar_file");

        if(scar_file_big) {
            print_to_stdout_and_file("Loading fibrosis patterns from file %s\n", scar_file_big);
            set_human_mesh_fibrosis_from_file(the_grid, 'b', scar_file_big, 2172089);
        }

        if(scar_file_small) {
            print_to_stdout_and_file("Loading fibrosis patterns from file %s\n", scar_file_small);
            set_human_mesh_fibrosis_from_file(the_grid, 's', scar_file_small, 845051);
        }

        if(!(scar_file_big || scar_file_small)) {

            real_cpu small_scar_center_x = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_x, config->config_data.config,
                                                        "small_scar_center_x");

            real_cpu small_scar_center_y = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_y, config->config_data.config,
                                                        "small_scar_center_y");

            real_cpu small_scar_center_z = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_z, config->config_data.config,
                                                        "small_scar_center_z");

            real_cpu big_scar_center_x = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_x, config->config_data.config,
                                                        "big_scar_center_x");

            real_cpu big_scar_center_y = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_y, config->config_data.config,
                                                        "big_scar_center_y");

            real_cpu big_scar_center_z = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_z, config->config_data.config,
                                                        "big_scar_center_z");

            real_cpu phi = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data.config, "phi");

            unsigned seed = 0;
            bool seed_success;
            GET_PARAMETER_NUMERIC_VALUE(unsigned, seed, config->config_data.config, "seed", seed_success);
            if(!seed_success)
                seed = 0;

            print_to_stdout_and_file("Setting random fibrosis pattern\n");
            set_human_mesh_fibrosis(the_grid, phi, seed, big_scar_center_x, big_scar_center_y, big_scar_center_z,
                                    small_scar_center_x, small_scar_center_y, small_scar_center_z);
        }
    }

    free(mesh_file);
    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_scar_wedge) {
    char *mesh_file;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(mesh_file, config->config_data.config, "mesh_file");

    char *scar_size;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(scar_size, config->config_data.config, "scar_size");

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data.config, "phi");

    unsigned fib_seed = 0;
    bool success;
    GET_PARAMETER_NUMERIC_VALUE(unsigned, fib_seed, config->config_data.config, "seed", success);

    if(!success)
        fib_seed = (unsigned)time(NULL) + getpid();

    srand(fib_seed);

    config->start_dx = 800.0;
    config->start_dy = 800.0;
    config->start_dz = 800.0;
    uint8_t size_code;

    initialize_and_construct_grid(the_grid, 204800, 204800, 204800);
    refine_grid(the_grid, 7);

    if(strcmp(scar_size, "big") == 0) {
        print_to_stdout_and_file("Loading Human Heart Edge with big scar\n");
        set_custom_mesh_with_bounds(the_grid, mesh_file, 2025252, 79100, 121000, 66700, 106000, 11200, 61400, true);
        size_code = 0;
    } else if(strcmp(scar_size, "small") == 0) {
        print_to_stdout_and_file("Loading Human Heart Edge with small scar\n");
        set_custom_mesh_with_bounds(the_grid, mesh_file, 2025252, 30400, 81600, 59200, 103000, 13600, 48000, true);
        size_code = 1;
    } else {
        print_to_stderr_and_file_and_exit("Function: initialize_grid_with_scar_edge, invalid scar size %s. Valid sizes are big or small. Exiting!\n",
            scar_size);
    }

    print_to_stdout_and_file("Cleaning grid\n");
    int i;
    for(i = 0; i < 7; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    refine_fibrotic_cells(the_grid);
    refine_fibrotic_cells(the_grid);
    refine_fibrotic_cells(the_grid);

    refine_border_zone_cells(the_grid);
    refine_border_zone_cells(the_grid);
    refine_border_zone_cells(the_grid);

    real_cpu scar_center_x;
    real_cpu scar_center_y;
    real_cpu scar_center_z;

    ////Fibrosis configuration

    // BIG SCAR
    if(size_code == 0) {
        scar_center_x = 95300.0;
        scar_center_y = 81600.0;
        scar_center_z = 36800.0;
    } else {
        scar_center_x = 52469.0;
        scar_center_y = 83225.0;
        scar_center_z = 24791.0;
    }

    real_cpu bz_size = 0.0;
    real_cpu dist;

    print_to_stdout_and_file("Using %u as seed\n", fib_seed);
    print_to_stdout_and_file("Calculating fibrosis using phi: %lf\n", phi);
    struct cell_node *grid_cell = the_grid->first_cell;
    bool fibrotic, border_zone;

    while(grid_cell != 0) {

        if(grid_cell->active) {
            fibrotic = FIBROTIC(grid_cell);
            border_zone = BORDER_ZONE(grid_cell);

            if(fibrotic) {
                grid_cell->can_change = false;
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX); // rand() has limited randomness
                if(p < phi)
                    grid_cell->active = false;
            } else if(border_zone) {
                real_cpu center_x = grid_cell->center_x;
                real_cpu center_y = grid_cell->center_y;
                real_cpu center_z = grid_cell->center_z;
                dist = sqrt((center_x - scar_center_x) * (center_x - scar_center_x) +
                            (center_y - scar_center_y) * (center_y - scar_center_y) +
                            (center_z - scar_center_z) * (center_z - scar_center_z));
                if(dist > bz_size) {
                    bz_size = dist;
                }
            }
        }
        grid_cell = grid_cell->next;
    }

    grid_cell = the_grid->first_cell;
    while(grid_cell != 0) {

        if(grid_cell->active) {
            border_zone = BORDER_ZONE(grid_cell);
            if(border_zone) {
                real_cpu center_x = grid_cell->center_x;
                real_cpu center_y = grid_cell->center_y;
                real_cpu center_z = grid_cell->center_z;
                dist = sqrt((center_x - scar_center_x) * (center_x - scar_center_x) +
                            (center_y - scar_center_y) * (center_y - scar_center_y) +
                            (center_z - scar_center_z) * (center_z - scar_center_z));
                dist = dist / bz_size;

                real_cpu phi_local = phi - phi * dist;
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);

                if(p < phi_local)
                    grid_cell->active = false;

                grid_cell->can_change = false;
            }
        }
        grid_cell = grid_cell->next;
    }

    free(mesh_file);
    free(scar_size);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_rabbit_mesh) {

    config->start_dx = 250.0;
    config->start_dy = 250.0;
    config->start_dz = 250.0;

    char *mesh_file = strdup("meshes/rabheart.alg");
    //LEAK on mesh_file if mesh_file is defined on ini file
    GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(mesh_file, config->config_data.config, "mesh_file");

    real_cpu max_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_h, config->config_data.config, "maximum_discretization");

    initialize_and_construct_grid(the_grid, 64000.0, 64000.0, 64000.0);
    refine_grid(the_grid, 7);

    print_to_stdout_and_file("Loading Rabbit Heart Mesh\n");

    set_custom_mesh(the_grid, mesh_file, 470197, false);

    print_to_stdout_and_file("Cleaning grid\n");
    int i;
    for(i = 0; i < 6; i++) {
        derefine_grid_inactive_cells(the_grid);
    }
    free(mesh_file);

    config->max_dx = max_h;
    config->max_dx_was_set = true;

    config->max_dy = max_h;
    config->max_dy_was_set = true;

    config->max_dz = max_h;
    config->max_dz_was_set = true;

    config->start_dx = 250.0;
    config->start_dx_was_set = true;
    config->start_dy = 250.0;
    config->start_dy_was_set = true;
    config->start_dz = 250.0;
    config->start_dz_was_set = true;

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_mouse_mesh) {

    char *mesh_file = NULL;

    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(mesh_file, config->config_data.config, "mesh_file");


    real_cpu start_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_h, config->config_data.config, "start_discretization");

    real_cpu max_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_h, config->config_data.config, "maximum_discretization");

    assert(the_grid);

    initialize_and_construct_grid(the_grid, 6400.0, 6400.0, 6400.0);

    refine_grid(the_grid, 5);

    print_to_stdout_and_file("Loading Mouse Heart Mesh\n");

    set_custom_mesh(the_grid, mesh_file, 96195, false);

    int i;
    for(i = 0; i < 5; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    if(start_h == 100.0) {

    } else if(start_h == 50.0) {
        print_to_stdout_and_file("Refining Mesh to 50um\n");
        refine_grid(the_grid, 1);
    } else if(start_h == 25.0) {
        print_to_stdout_and_file("Refining Mesh to 25um\n");
        refine_grid(the_grid, 2);
    } else if(start_h == 12.5) {
        print_to_stdout_and_file("Refining Mesh to 12.5um\n");
        refine_grid(the_grid, 3);
    } else {
        print_to_stderr_and_file_and_exit("Invalid discretizations for this mesh. Valid discretizations are: 100um, 50um, 25um "
                                 "or 12.5um. Using 100um!\n");
        start_h = 100.0;
    }

    config->max_dx = max_h;
    config->max_dx_was_set = true;

    config->max_dy = max_h;
    config->max_dy_was_set = true;

    config->max_dz = max_h;
    config->max_dz_was_set = true;

    config->start_dx = start_h;
    config->start_dx_was_set = true;
    config->start_dy = start_h;
    config->start_dy_was_set = true;
    config->start_dz = start_h;
    config->start_dz_was_set = true;

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_benchmark_mesh) {

    real_cpu side_length;

    real_cpu start_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_h, config->config_data.config, "start_discretization");

    real_cpu max_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_h, config->config_data.config, "maximum_discretization");

    print_to_stdout_and_file("Loading N-Version benchmark mesh using dx %lf um, dy %lf um, dz %lf um\n", start_h,
                             start_h, start_h);

    side_length = start_h;

    while(side_length < 20000.0) {
        side_length = side_length * 2.0;
    }

    initialize_and_construct_grid(the_grid, side_length, side_length, side_length);

    config->max_dx = max_h;
    config->max_dx_was_set = true;

    config->max_dy = max_h;
    config->max_dy_was_set = true;

    config->max_dz = max_h;
    config->max_dz_was_set = true;

    config->start_dx = start_h;
    config->start_dx_was_set = true;
    config->start_dy = start_h;
    config->start_dy_was_set = true;
    config->start_dz = start_h;
    config->start_dz_was_set = true;


    int num_steps = get_num_refinement_steps_to_discretization(side_length, start_h);

    refine_grid(the_grid, num_steps);
    set_benchmark_domain(the_grid);

    print_to_stdout_and_file("Cleaning grid\n");
    int i;

    for(i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    if(the_grid->adaptive) {
        struct cell_node *grid_cell;
        grid_cell = the_grid->first_cell;

        while(grid_cell != 0) {
            if(grid_cell->active) {
                set_cell_not_changeable(grid_cell, start_h);
            }
            grid_cell = grid_cell->next;
        }
    }

    translate_visible_mesh_to_origin(the_grid);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_fibrotic_mesh) {

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data.config, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config->config_data.config, "seed");

    initialize_grid_with_square_mesh(config, the_grid);
    set_plain_fibrosis(the_grid, phi, seed);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_source_sink_fibrotic_mesh) 
{

    real_cpu channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, channel_width, config->config_data.config, "channel_width");

    real_cpu channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, channel_length, config->config_data.config, "channel_length");


    initialize_grid_with_square_mesh(config, the_grid);
    set_plain_source_sink_fibrosis(the_grid, channel_width, channel_length);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_and_sphere_fibrotic_mesh) {

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data.config, "phi");

    real_cpu plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, plain_center, config->config_data.config, "plain_center");

    real_cpu sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sphere_radius, config->config_data.config, "sphere_radius");

    real_cpu border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, border_zone_size, config->config_data.config,
                                                "border_zone_size");

    real_cpu border_zone_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, border_zone_radius, config->config_data.config,
                                                "border_zone_radius");

    bool success;
    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE(unsigned, seed, config->config_data.config, "seed", success);
    if(!success) {
        seed = 0;
    }

    initialize_grid_with_square_mesh(config, the_grid);
    set_plain_sphere_fibrosis(the_grid, phi, plain_center, sphere_radius, border_zone_size, border_zone_radius, seed);

    return 1;
}
