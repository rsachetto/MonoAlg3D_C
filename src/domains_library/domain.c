//
// Created by sachetto on 01/10/17.
//

#include "domain_helpers.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"
#include "../utils/utils.h"
#include <time.h>
#include <unistd.h>

SET_SPATIAL_DOMAIN(initialize_grid_with_cuboid_mesh) {

    real_cpu start_dx = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dx, config, "start_dx");

    real_cpu start_dy = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dy, config, "start_dy");

    real_cpu start_dz = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dz, config, "start_dz");

    real_cpu side_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length_x, config, "side_length_x");

    real_cpu side_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length_y, config, "side_length_y");

    real_cpu side_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length_z, config, "side_length_z");

    return set_cuboid_domain_mesh(the_grid, start_dx, start_dy, start_dz, side_length_x, side_length_y, side_length_z);

}

SET_SPATIAL_DOMAIN(initialize_grid_with_spherical_mesh) {

    real_cpu diameter = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, diameter, config, "diameter");

    real_cpu start_dx = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dx, config, "start_dx");

    real_cpu start_dy = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dy, config, "start_dy");

    real_cpu start_dz = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dz, config, "start_dz");

    real_cpu radius = diameter / 2.0;

    int ret = set_cuboid_domain_mesh(the_grid, start_dx, start_dy, start_dz, diameter, diameter, diameter);

    if(ret == 0) {
        return 0;
    }

    // Distance between the sphere center and a cell.
    double distance;

    // Coordinates of the sphere center
    const real_cpu x_c = 1050.0;
    const real_cpu y_c = 1050.0;
    const real_cpu z_c = 1050.0;

    double x, y, z;
    FOR_EACH_CELL(the_grid) {
        x = cell->center.x;
        y = cell->center.y;
        z = cell->center.z;
        distance = (x - x_c) * (x - x_c) + (y - y_c) * (y - y_c) + (z - z_c) * (z - z_c);

        cell->active = (distance <= (radius * radius));
    }

    return ret;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_square_mesh) {
    return set_square_mesh(config, the_grid);
}

SET_SPATIAL_DOMAIN(initialize_grid_with_cable_mesh) {

    real_cpu start_dx = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dx, config, "start_dx");

    real_cpu start_dy = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dy, config, "start_dy");

    real_cpu start_dz = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dz, config, "start_dz");

    real_cpu cable_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, cable_length, config, "cable_length");

    return set_cuboid_domain_mesh(the_grid, start_dx, start_dy, start_dz, cable_length, start_dy, start_dz);

}

SET_SPATIAL_DOMAIN(initialize_grid_with_rabbit_mesh) {

    real_cpu start_discretization = 250;

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    real_cpu maximum_discretization = start_discretization;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, maximum_discretization, config, "maximum_discretization");

    log_info("Loading Rabbit Heart Mesh\n");

    uint32_t num_volumes = 470197;
    uint32_t num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_discretization, 0, NULL);

    log_info("Read %d volumes from file: %s\n", num_loaded, mesh_file);

    free(mesh_file);

    the_grid->start_discretization = SAME_POINT3D(start_discretization);
    the_grid->max_discretization = SAME_POINT3D(maximum_discretization);

    return num_loaded > 0;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_benchmark_mesh) {

    real_cpu side_length;

    real_cpu start_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_h, config, "start_discretization");

    real_cpu max_h = start_h;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_h, config, "maximum_discretization");

    log_info("Loading N-Version benchmark mesh using dx %lf um, dy %lf um, dz %lf um\n", start_h, start_h, start_h);

    side_length = start_h;

    while(side_length < 20000.0) {
        side_length = side_length * 2.0;
    }

    initialize_and_construct_grid(the_grid, POINT3D(side_length, side_length, side_length));

    int num_steps = get_num_refinement_steps_to_discretization(side_length, start_h);

    refine_grid(the_grid, num_steps);
    set_benchmark_domain(the_grid);

    log_info("Cleaning grid\n");
    int i;

    for(i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    if(the_grid->adaptive) {
        FOR_EACH_CELL(the_grid) {
            if(cell->active) {
                set_cell_not_changeable(cell, start_h);
            }
        }
    }

    the_grid->start_discretization = SAME_POINT3D(start_h);
    the_grid->max_discretization = SAME_POINT3D(max_h);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_fibrotic_mesh) {

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config, "seed");

    set_square_mesh(config, the_grid);

    if(seed == 0)
        seed = (unsigned)time(NULL) + getpid();

    srand(seed);

    ADD_UINT_PARAMETER_TO_CONFIG("seed", seed, config);

    set_plain_fibrosis(the_grid, phi, seed);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_fibrotic_mesh_from_file) {

    char *fib_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(fib_file, config, "fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config, "size");

    set_square_mesh(config, the_grid);
    set_fibrosis_from_file(the_grid, fib_file, fib_size);

    free(fib_file);
    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_source_sink_fibrotic_mesh) {

    real_cpu channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, channel_width, config, "channel_width");

    real_cpu channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, channel_length, config, "channel_length");

    set_square_mesh(config, the_grid);
    set_plain_source_sink_fibrosis(the_grid, channel_width, channel_length);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_and_sphere_fibrotic_mesh) {

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config, "phi");

    real_cpu plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, plain_center, config, "plain_center");

    real_cpu sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sphere_radius, config, "sphere_radius");

    real_cpu border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, border_zone_size, config, "border_zone_size");

    real_cpu border_zone_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, border_zone_radius, config, "border_zone_radius");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config, "seed");

    set_square_mesh(config, the_grid);
    set_plain_sphere_fibrosis(the_grid, phi, plain_center, sphere_radius, border_zone_size, border_zone_radius, seed);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_cuboid_and_sphere_fibrotic_mesh) {

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config, "phi");

    real_cpu sphere_center[3] = {0,0,0};
    GET_PARAMETER_VECTOR3_VALUE_OR_USE_DEFAULT(sphere_center, config, "sphere_center");

    real_cpu sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sphere_radius, config, "sphere_radius");


    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config, "seed");

    initialize_grid_with_cuboid_mesh(config, the_grid);
    set_cube_sphere_fibrosis(the_grid, phi, sphere_center, sphere_radius, seed);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_and_sphere_fibrotic_mesh_without_inactivating) {

    real_cpu plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, plain_center, config, "plain_center");

    real_cpu sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sphere_radius, config, "sphere_radius");

    real_cpu border_zone_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, border_zone_radius, config, "border_zone_radius");

    set_square_mesh(config, the_grid);

    set_plain_sphere_fibrosis_without_inactivating(the_grid, plain_center, sphere_radius, border_zone_radius);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_square_mesh_and_fibrotic_region) {

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config, "seed");

    real min_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_x, config, "region_min_x");

    real max_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_x, config, "region_max_x");

    real min_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_y, config, "region_min_y");

    real max_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_y, config, "region_max_y");

    real min_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_z, config, "region_min_z");

    real max_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_z, config, "region_max_z");

    set_square_mesh(config, the_grid);
    set_plain_fibrosis_inside_region(the_grid, phi, seed, min_x, max_x, min_y, max_y, min_z, max_z);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_custom_mesh) {

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    real_cpu start_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_h, config, "start_discretization");

    real_cpu max_h = start_h;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_h, config, "maximum_discretization");

    uint32_t total_number_mesh_points = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, total_number_mesh_points, config, "number_of_points");

    the_grid->start_discretization = SAME_POINT3D(start_h);
    the_grid->max_discretization = SAME_POINT3D(max_h);

    int ret = (int) set_custom_mesh_from_file(the_grid, mesh_file, total_number_mesh_points, start_h, 0, NULL);

    free(mesh_file);

    return ret;

}
