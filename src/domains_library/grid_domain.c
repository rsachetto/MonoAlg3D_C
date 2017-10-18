//
// Created by sachetto on 01/10/17.
//

#include "grid_domain.h"

#include "../hash/point_hash.h"
#include "../utils/utils.h"
#include "../utils/erros_helpers.h"
#include "../utils/logfile_utils.h"
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>


int get_num_refinement_steps_to_discretization (double side_len, double h) {

    int num_steps = 0;
    double aux = side_len;

    while (aux > h) {
        num_steps++;
        aux /= 2.0;
    }

    return num_steps - 1;
}

/**
 * Sets the current domain as a domain described in the N-version benchmark
 * (http://rsta.royalsocietypublishing.org/content/369/1954/4331)
 *
 */
void set_benchmark_domain (struct grid *the_grid) {
    struct cell_node *grid_cell = the_grid->first_cell;
    while (grid_cell != 0) {
        grid_cell->active = (grid_cell->center_y < 20000) && (grid_cell->center_x < 7000) &&
                            (grid_cell->center_z < 3000);
        grid_cell = grid_cell->next;
    }
}

void set_plain_domain (struct grid *the_grid, double sizeX, double sizeY, double sizeZ) {
    struct cell_node *grid_cell = the_grid->first_cell;

    while (grid_cell != 0) {
        grid_cell->active = (grid_cell->center_y < sizeY) && (grid_cell->center_x < sizeX) &&
                            (grid_cell->center_z < sizeZ);
        grid_cell = grid_cell->next;
    }
}

void set_custom_mesh (struct grid *the_grid, const char *file_name, int size) {
    struct cell_node *grid_cell = the_grid->first_cell;
    FILE *file = fopen (file_name, "r");

    if (!file) {
        print_to_stdout_and_file ("Error opening mesh described in %s!!\n", file_name);
        exit (0);
    }

    double **a = (double **)malloc (sizeof (double *) * size);
    for (int i = 0; i < size; i++) {
        a[i] = (double *)malloc (sizeof (double) * 3);
        if (a[i] == NULL) {
            print_to_stdout_and_file ("Failed to allocate memory\n");
            exit (0);
        }
    }
    double b;

    double maxy = 0.0;
    double maxz = 0.0;
    double miny = DBL_MAX;
    double minz = DBL_MAX;
    int fibrosis;

    int i = 0;

    struct point_hash *p_hash = point_hash_create ();
    struct point_3d point3d;

    while (i < size) {
        fscanf (file, "%lf,%lf,%lf,%lf,%d\n", &a[i][0], &a[i][1], &a[i][2], &b, &fibrosis);
        if (a[i][1] > maxy)
            maxy = a[i][1];
        if (a[i][2] > maxz)
            maxz = a[i][2];
        if (a[i][1] < miny)
            miny = a[i][1];
        if (a[i][2] < minz)
            minz = a[i][2];

        point3d.x = a[i][0];
        point3d.y = a[i][1];
        point3d.z = a[i][2];

        point_hash_insert (p_hash, point3d, fibrosis);

        i++;
    }
    sort_vector (a, size);

    double maxx = a[size - 1][0];
    double minx = a[0][0];

    double x, y, z;
    while (grid_cell != 0) {
        x = grid_cell->center_x;
        y = grid_cell->center_y;
        z = grid_cell->center_z;

        point3d.x = x;
        point3d.y = y;
        point3d.z = z;

        grid_cell->fibrotic = (point_hash_search (p_hash, point3d) == 1);
        grid_cell->border_zone = (point_hash_search (p_hash, point3d) == 2);

        if (x > maxx || y > maxy || z > maxz || x < minx || y < miny || z < minz) {

            grid_cell->active = false;
        } else {
            int index = inside_mesh (a, x, y, z, 0, size - 1);
            grid_cell->active = index != -1;
        }
        grid_cell = grid_cell->next;
    }
    fclose (file);
    // deallocate memory
    int l;
    for (l = 0; l < size; l++) {
        free (a[l]);
    }
    free (a);
    point_hash_destroy (p_hash);
}

void set_custom_mesh_with_bounds (struct grid *the_grid, const char *file_name, int size,
                                  double minx, double maxx, double miny, double maxy, double minz,
                                  double maxz) {
    struct cell_node *grid_cell = the_grid->first_cell;
    FILE *file = fopen (file_name, "r");

    if (!file) {
        print_to_stdout_and_file ("Error opening mesh described in %s!!\n", file_name);
        exit (0);
    }

    double **a = (double **)malloc (sizeof (double *) * size);
    for (int i = 0; i < size; i++) {
        a[i] = (double *)malloc (sizeof (double) * 3);
        if (a[i] == NULL) {
            print_to_stdout_and_file ("Failed to allocate memory\n");
            exit (0);
        }
    }
    double b;

    int fibrosis;

    int i = 0;

    struct point_hash *p_hash = point_hash_create ();
    struct point_3d point3d;

    while (i < size) {
        fscanf (file, "%lf,%lf,%lf,%lf,%d\n", &a[i][0], &a[i][1], &a[i][2], &b, &fibrosis);

        point3d.x = a[i][0];
        point3d.y = a[i][1];
        point3d.z = a[i][2];

        // myMap[TPoint3D(a[i])] = fibrosis;
        point_hash_insert (p_hash, point3d, fibrosis);

        i++;
    }
    sort_vector (a, size);

    double x, y, z;
    while (grid_cell != 0) {
        x = grid_cell->center_x;
        y = grid_cell->center_y;
        z = grid_cell->center_z;

        point3d.x = x;
        point3d.y = y;
        point3d.z = z;

        grid_cell->fibrotic = (point_hash_search (p_hash, point3d) == 1);
        grid_cell->border_zone = (point_hash_search (p_hash, point3d) == 2);

        if (x > maxx || y > maxy || z > maxz || x < minx || y < miny || z < minz) {

            grid_cell->active = false;
        } else {
            int index = inside_mesh (a, x, y, z, 0, size - 1);
            grid_cell->active = index != -1;
        }
        grid_cell = grid_cell->next;
    }
    fclose (file);
    // deallocate memory
    int j;
    for (j = 0; j < size; j++) {
        free (a[j]);
    }
    free (a);
    point_hash_destroy (p_hash);
}

void set_human_mesh (struct grid *the_grid, const char *file_name) {
    set_custom_mesh (the_grid, file_name, 2025252);
}

void set_human_sub_mesh (struct grid *the_grid, const char *file_name, double minx, double maxx,
                         double miny, double maxy, double minz, double maxz) {
    set_custom_mesh_with_bounds (the_grid, file_name, 2025252, minx, maxx, miny, maxy, minz, maxz);
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

    set_custom_mesh (the_grid, mesh_file, 470197);

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
    
    set_custom_mesh (the_grid, mesh_file, 96195);

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

void set_cell_not_changeable(struct cell_node *c, double initialDiscretization) {

    double P1x, P1y, P1z;
    double P2x, P2y, P2z;
    double P3x, P3y, P3z;
    double P4x, P4y, P4z;
    double P5x, P5y, P5z;
    double P6x, P6y, P6z;
    double P7x, P7y, P7z;
    double P8x, P8y, P8z;
    double Cx, Cy, Cz;

    if(initialDiscretization == 100.0) {
        P1x = 6950; P1y = 50;    P1z = 50;
        P2x = 6950; P2y = 19950; P2z = 50;
        P3x = 6950; P3y = 50;    P3z = 2950;
        P4x = 6950; P4y = 19950; P4z = 2950;
        P5x = 50;   P5y = 50 ;   P5z = 50;
        P6x = 50;   P6y = 19950; P6z = 50;
        P7x = 50;   P7y = 50;    P7z = 2950;
        P8x = 50;   P8y = 19950; P8z = 2950;
        Cx = 3450;  Cy  = 9950;    Cz = 1450;
    }

    else if (initialDiscretization == 200.0) {
        P1x = 6900; P1y = 100; P1z = 100;
        P2x = 6900; P2y = 19900; P2z = 100;
        P3x = 6900; P3y = 100; P3z = 2900;
        P4x = 6900; P4y = 19900; P4z = 2900;
        P5x = 100; P5y = 100; P5z = 100;
        P6x = 100; P6y = 19900; P6z = 100;
        P7x = 100; P7y = 100; P7z = 2900;
        P8x = 100; P8y = 19900; P8z = 2900;
        Cx = 3500; Cy = 9900; Cz = 1500;
    }

    else if (initialDiscretization == 125.0) {
        P1x = 6937.5; P1y = 62.5; P1z = 62.5;
        P2x = 6937.5; P2y = 19937.5; P2z = 62.5;
        P3x = 6937.5; P3y = 62.5; P3z = 2937.5;
        P4x = 6937.5; P4y = 19937.5; P4z = 2937.5;
        P5x = 62.5; P5y = 62.5; P5z = 62.5;
        P6x = 62.5; P6y = 19937.5; P6z = 62.5;
        P7x = 3937.5; P7y = 19937.5; P7z = 62.5;
        P8x = 62.5; P8y = 19937.5; P8z = 2937.5;
        Cx = 3437.5; Cy = 9937.5; Cz = 1562.5;
    }

    else if (initialDiscretization == 250.0) {
        P1x = 6875; P1y = 125; P1z = 125;
        P2x = 6875; P2y = 19875; P2z = 125;
        P3x = 6875; P3y = 125; P3z = 2875;
        P4x = 6875; P4y = 19875; P4z = 2875;
        P5x = 125; P5y = 125; P5z = 125;
        P6x = 125; P6y = 19875; P6z = 125;
        P7x = 125; P7y = 125; P7z = 2875;
        P8x = 125; P8y = 19875; P8z = 2875;
        Cx = 3375; Cy = 9875; Cz = 1125;
    }

    else {
        P1x = -1; P1y = -1; P1z = -1;
        P2x = -1; P2y = -1; P2z = -1;
        P3x = -1; P3y = -1; P3z = -1;
        P4x = -1; P4y = -1; P4z = -1;
        P5x = -1; P5y = -1; P5z = -1;
        P6x = -1; P6y = -1; P6z = -1;
        P7x = -1; P7y = -1; P7z = -1;
        P8x = -1; P8y = -1; P8z = -1;
        Cx = -1; Cy = -1; Cz = -1;
    }
    bool cannotChange = ( ( c->center_x == P1x ) && ( c->center_y == P1y ) && ( c->center_z == P1z ) );
    cannotChange |= ( ( c->center_x == P2x ) && ( c->center_y == P2y ) && ( c->center_z == P2z ) );
    cannotChange |= ( ( c->center_x == P3x ) && ( c->center_y == P3y ) && ( c->center_z == P3z ) );
    cannotChange |= ( ( c->center_x == P4x ) && ( c->center_y == P4y ) && ( c->center_z == P4z ) );
    cannotChange |= ( ( c->center_x == P5x ) && ( c->center_y == P5y ) && ( c->center_z == P5z ) );
    cannotChange |= ( ( c->center_x == P6x ) && ( c->center_y == P6y ) && ( c->center_z == P6z ) );
    cannotChange |= ( ( c->center_x == P7x ) && ( c->center_y == P7y ) && ( c->center_z == P7z ) );
    cannotChange |= ( ( c->center_x == P8x ) && ( c->center_y == P8y ) && ( c->center_z == P8z ) );
    cannotChange |= ( ( c->center_x == Cx )  && ( c->center_y == Cy )  && ( c->center_z == Cz ) );

    c->can_change = !cannotChange;



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

void set_plain_fibrosis(struct grid* the_grid, double phi, unsigned fib_seed) {

    print_to_stdout_and_file("Making %.2lf %% of cells inactive\n", phi*100.0);

    struct cell_node *grid_cell;

    if(fib_seed == 0)
        fib_seed = (unsigned) time(NULL) + getpid();

    srand(fib_seed);


    print_to_stdout_and_file("Using %u as seed\n", fib_seed);

    grid_cell = the_grid->first_cell;
    while( grid_cell != 0 ) {

        if(grid_cell->active) {
            double p = (double) (rand()) / (RAND_MAX);
            if (p < phi) grid_cell->active = false;
            grid_cell->fibrotic = true;

        }
        grid_cell = grid_cell->next;
    }

}

void set_plain_sphere_fibrosis(struct grid* the_grid, double phi,  double plain_center, double sphere_radius, double bz_size,
                               double bz_radius,  unsigned fib_seed) {

    print_to_stdout_and_file("Making %.2lf %% of cells inactive\n", phi*100.0f);

    if(fib_seed == 0)
        fib_seed = (unsigned) time(NULL) + getpid();

    srand(fib_seed);

    print_to_stdout_and_file("Using %u as seed\n", fib_seed);

    double bz_radius_2 = pow(bz_radius, 2.0);
    double sphere_radius_2 = pow(sphere_radius, 2.0);

    struct cell_node *grid_cell;

    grid_cell = the_grid->first_cell;
    while (grid_cell != 0) {

        double distance = pow(grid_cell->center_x - plain_center, 2.0) +
                          pow(grid_cell->center_y - plain_center, 2.0);

        if (grid_cell->active) {
            if (distance <= bz_radius_2) {
                if (distance <= sphere_radius_2) {
                    grid_cell->fibrotic = true;
                } else {
                    grid_cell->border_zone = true;
                }
            }
        }
        grid_cell = grid_cell->next;
    }

    grid_cell = the_grid->first_cell;

    while (grid_cell != 0) {

        if (grid_cell->active) {
            if (grid_cell->fibrotic) {
                double p = (double) (rand()) / (RAND_MAX);
                if (p < phi)
                    grid_cell->active = false;
                grid_cell->can_change = false;
            } else if (grid_cell->border_zone) {
                double distance_from_center = sqrt(
                        (grid_cell->center_x - plain_center) * (grid_cell->center_x - plain_center) +
                        (grid_cell->center_y - plain_center) * (grid_cell->center_y - plain_center));
                distance_from_center = (distance_from_center - sphere_radius) / bz_size;
                double phi_local = phi - phi * distance_from_center;
                double p = (double) (rand()) / (RAND_MAX);
                if (p < phi_local)
                    grid_cell->active = false;
                grid_cell->can_change = false;
            }
        }
        grid_cell = grid_cell->next;
    }


}