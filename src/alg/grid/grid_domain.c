//
// Created by sachetto on 01/10/17.
//

#include "../../utils/hash/point_hash.h"
#include "../../utils/utils.h"
#include "grid.h"
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

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
        printf ("Error opening mesh described in %s!!\n", file_name);
        exit (0);
    }

    double **a = (double **)malloc (sizeof (double *) * size);
    for (int i = 0; i < size; i++) {
        a[i] = (double *)malloc (sizeof (double) * 3);
        if (a[i] == NULL) {
            printf ("Failed to allocate memory\n");
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
        printf ("Error opening mesh described in %s!!\n", file_name);
        exit (0);
    }

    double **a = (double **)malloc (sizeof (double *) * size);
    for (int i = 0; i < size; i++) {
        a[i] = (double *)malloc (sizeof (double) * 3);
        if (a[i] == NULL) {
            printf ("Failed to allocate memory\n");
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

void initialize_grid_with_rabbit_mesh (struct grid *the_grid, const char *mesh_file) {

    initialize_and_construct_grid (the_grid, 32000.0, the_grid->num_cell_neighbours);
    refine_grid (the_grid, 6);

    printf ("Loading Rabbit Heart Mesh\n");
    ;
    set_custom_mesh (the_grid, mesh_file, 470197);

    printf ("Cleaning grid\n");
    int i;
    for (i = 0; i < 6; i++) {
        derefine_grid_inactive_cells (the_grid);
    }
}

void initialize_grid_with_mouse_mesh (struct grid *the_grid, const char *mesh_file) {

    if (the_grid == NULL) {
        fprintf (stderr, "set_mouse_mesh: the_grid is NULL. Exiting!");
        exit (10);
    }

    initialize_and_construct_grid (the_grid, 6400.0, 7);

    refine_grid (the_grid, 5);

    printf ("Loading Mouse Heart Mesh\n");
    ;
    set_custom_mesh (the_grid, mesh_file, 96195);

    printf ("Cleaning grid\n");

    int i;
    for (i = 0; i < 5; i++) {
        derefine_grid_inactive_cells (the_grid);
    }
}

void initialize_grid_with_benchmark_mesh (struct grid *the_grid, double start_h) {

    double side_length;

    printf ("Loading N-Version benchmark mesh using dx %lf um\n", start_h);
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

    printf ("Cleaning grid\n");
    int i;

    for (i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells (the_grid);
    }

    //if (the_grid->adaptive) {
        // TODO: @incomplete: set fixed cells for activation time caculation (maybe a separate
        // funcition??)
        //        cout << "Setting fixed cells for Activation time calculation" << endl;
        //        CellNode *grid_cell;
        //        grid_cell = firstCell;
        //
        //        while( grid_cell != 0 ) {
        //            if(grid_cell->active) {
        //                setCellNotChangeable(grid_cell, globalArgs.start_h);
        //            }
        //            grid_cell = grid_cell->next;
        //        }
    //}
    // exit(0);
}

void initialize_grid_with_plain_mesh (struct grid *the_grid, double desired_side_lenght, double start_h, int num_layers) {

    double real_side_lenght = start_h * 2.0f;
    double max_h = start_h * num_layers;

    while (real_side_lenght < desired_side_lenght) {
        real_side_lenght *= 2.0f;
    }

    printf("Initial cube side length: %lfum x %lfum x %lfum\n", real_side_lenght, real_side_lenght, real_side_lenght);
    printf("Loading plain mesh with %lfum x %lfum x %lfum using dx %lfum\n", desired_side_lenght, desired_side_lenght, max_h, start_h);

    int num_steps = get_num_refinement_steps_to_discretization(real_side_lenght, start_h);

    initialize_and_construct_grid(the_grid, real_side_lenght, 7);

    if ((real_side_lenght / 2.0f) > max_h) {
        double aux = real_side_lenght/2.0f;

        for (int i = 0; i < num_steps - 3; i++) {
            set_plain_domain(the_grid, real_side_lenght, real_side_lenght, aux);
            refine_grid(the_grid, 1);
            aux = aux / 2.0f;
        }

        refine_grid(the_grid, 3);


    } else {
        refine_grid(the_grid, num_steps);
    }

    set_plain_domain(the_grid, desired_side_lenght, desired_side_lenght, max_h);

    int i;
    for (i = 0; i < num_steps; i++) {
        derefine_grid_inactive_cells(the_grid);
    }


}

void initialize_grid_with_plain_fibrotic_mesh(struct grid *the_grid, double side_length, double start_h, int num_layers,
                                              double phi) {


    initialize_grid_with_plain_mesh(the_grid, side_length, start_h, num_layers);
    set_plain_fibrosis(the_grid, phi);



}


void initialize_grid_with_plain_and_sphere_fibrotic_mesh(struct grid *the_grid, double side_length,
                                                         double start_h, int num_layers, double phi,
                                                         double plain_center, double sphere_radius, double bz_size,
                                                         double bz_radius) {

    initialize_grid_with_plain_mesh(the_grid, side_length, start_h, num_layers);
    set_plain_sphere_fibrosis(the_grid, phi, plain_center,sphere_radius,bz_size, bz_radius);


}

void set_plain_fibrosis(struct grid* the_grid, double phi) {

    printf("Making %.2lf %% of cells inactive\n", phi*100.0);

    struct cell_node *grid_cell;

    unsigned fib_seed = (unsigned) time(NULL) + getpid();
    srand(fib_seed);

    printf("Using %u as seed\n", fib_seed);

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
                               double bz_radius) {

    printf("Making %.2lf %% of cells inactive\n", phi*100.0f);

    unsigned fib_seed = (unsigned) time(NULL) + getpid();
    srand(fib_seed);

    printf("Using %u as seed\n", fib_seed);

    double bz_radius_2 = pow(bz_radius, 2.0f);
    double sphere_radius_2 = pow(sphere_radius, 2.0f);

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

void save_grid_domain (struct grid * the_grid, const char *file_name) {
    struct cell_node *grid_cell = the_grid->first_cell;
    FILE *f = fopen (file_name, "w");

    while (grid_cell != 0) {
        if (grid_cell->active) {
            fprintf (f, "%lf,%lf,%lf,%lf\n", grid_cell->center_x, grid_cell->center_y,
                     grid_cell->center_z, grid_cell->half_face_length);
        }
        grid_cell = grid_cell->next;
    }
    fclose (f);
}
