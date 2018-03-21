//
// Created by sachetto on 19/10/17.
//

#include "domain_helpers.h"
#include "../utils/logfile_utils.h"
#include "../utils/utils.h"
#include <float.h>
#include <math.h>
#include <time.h>

#ifdef _MSC_VER
#include <process.h>
    #define getpid _getpid
#else
#include <unistd.h>
#endif

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
        grid_cell->active =
                (grid_cell->center_y < 20000) && (grid_cell->center_x < 7000) && (grid_cell->center_z < 3000);
        grid_cell = grid_cell->next;
    }
}

void set_plain_domain (struct grid *the_grid, double sizeX, double sizeY, double sizeZ) {
    struct cell_node *grid_cell = the_grid->first_cell;

    while (grid_cell != 0) {
        grid_cell->active =
                (grid_cell->center_y < sizeY) && (grid_cell->center_x < sizeX) && (grid_cell->center_z < sizeZ);
        grid_cell = grid_cell->next;
    }
}

void set_custom_mesh (struct grid *the_grid, const char *file_name, size_t size, bool read_fibrosis) {

    struct cell_node *grid_cell = the_grid->first_cell;
    FILE *file = fopen (file_name, "r");

    if (!file) {
        print_to_stdout_and_file ("Error opening mesh described in %s!!\n", file_name);
        exit (0);
    }

    double **mesh_points = (double **)malloc (sizeof (double *) * size);
    for (int i = 0; i < size; i++) {
        mesh_points[i] = (double *)malloc (sizeof (double) * 4);
        if (mesh_points[i] == NULL) {
            print_to_stdout_and_file ("Failed to allocate memory\n");
            exit (0);
        }
    }
    double dummy; // we don't use this value here

    double maxy = 0.0;
    double maxz = 0.0;
    double miny = DBL_MAX;
    double minz = DBL_MAX;
    int *fibrosis = (int *)malloc (sizeof (int) * size);

    char *tag = (char *)malloc (size);
    for (int k = 0; k < size; k++) {
        tag[k] = 'n';
    }

    int i = 0;
    while (i < size) {

        if(read_fibrosis) {
            fscanf (file, "%lf,%lf,%lf,%lf,%d,%c\n", &mesh_points[i][0], &mesh_points[i][1], &mesh_points[i][2], &dummy,
                    &fibrosis[i], &tag[i]);
        }
        else {
            fscanf (file, "%lf,%lf,%lf,%lf\n", &mesh_points[i][0], &mesh_points[i][1], &mesh_points[i][2], &dummy);
        }

        // we save the old index to reference fibrosis[i] and tags[i]. T
        // this is needed because the array mesh_points is sorted after reading the mesh file.
        mesh_points[i][3] = i;

        if (mesh_points[i][1] > maxy)
            maxy = mesh_points[i][1];
        if (mesh_points[i][2] > maxz)
            maxz = mesh_points[i][2];
        if (mesh_points[i][1] < miny)
            miny = mesh_points[i][1];
        if (mesh_points[i][2] < minz)
            minz = mesh_points[i][2];

        i++;
    }
    sort_vector (mesh_points, size); // we need to sort because inside_mesh perform a binary search

    double maxx = mesh_points[size - 1][0];
    double minx = mesh_points[0][0];
    int index;

    double x, y, z;
    while (grid_cell != 0) {
        x = grid_cell->center_x;
        y = grid_cell->center_y;
        z = grid_cell->center_z;

        if (x > maxx || y > maxy || z > maxz || x < minx || y < miny || z < minz) {
            grid_cell->active = false;
        } else {
            index = inside_mesh (mesh_points, x, y, z, 0, size - 1);

            if (index != -1) {
                grid_cell->active = true;
                if (read_fibrosis) {
                    int old_index = (int)mesh_points[index][3];
                    grid_cell->fibrotic = (fibrosis[old_index] == 1);
                    grid_cell->border_zone = (fibrosis[old_index] == 2);
                    grid_cell->scar_type = tag[old_index];
                }
            } else {
                grid_cell->active = false;
            }
        }
        grid_cell = grid_cell->next;
    }

    fclose (file);

    // deallocate memory
    for (int l = 0; l < size; l++) {
        free (mesh_points[l]);
    }

    free (mesh_points);
    free (tag);
    free (fibrosis);
}

void set_custom_mesh_with_bounds (struct grid *the_grid, const char *file_name, size_t size, double minx, double maxx,
                                  double miny, double maxy, double minz, double maxz, bool read_fibrosis) {

    struct cell_node *grid_cell = the_grid->first_cell;
    FILE *file = fopen (file_name, "r");

    if (!file) {
        print_to_stdout_and_file ("Error opening mesh described in %s!!\n", file_name);
        exit (0);
    }

    double **mesh_points = (double **)malloc (sizeof (double *) * size);
    for (int i = 0; i < size; i++) {
        mesh_points[i] = (double *)calloc (4, sizeof (double));
        if (mesh_points[i] == NULL) {
            print_to_stdout_and_file ("Failed to allocate memory\n");
            exit (0);
        }
    }
    double dummy; // we don't use this value here
    int *fibrosis = (int *)malloc (sizeof (int) * size);

    char *tag = (char *)malloc (size);
    for (int k = 0; k < size; k++) {
        tag[k] = 'n';
    }

    int i = 0;
    while (i < size) {

        fscanf (file, "%lf,%lf,%lf,%lf,%d,%c\n", &mesh_points[i][0], &mesh_points[i][1], &mesh_points[i][2], &dummy,
                &fibrosis[i], &tag[i]);

        // we save the old index to reference fibrosis[i] and tags[i]. T
        // this is needed because the array mesh_points is sorted after reading the mesh file.
        mesh_points[i][3] = i;
        i++;
    }
    sort_vector (mesh_points, size); // we need to sort because inside_mesh perform a binary search
    int index;

    double x, y, z;
    while (grid_cell != 0) {
        x = grid_cell->center_x;
        y = grid_cell->center_y;
        z = grid_cell->center_z;

        if (x > maxx || y > maxy || z > maxz || x < minx || y < miny || z < minz) {
            grid_cell->active = false;
        } else {
            index = inside_mesh (mesh_points, x, y, z, 0, size - 1);

            if (index != -1) {
                grid_cell->active = true;
                if (read_fibrosis) {
                    int old_index = (int)mesh_points[index][3];
                    grid_cell->fibrotic = (fibrosis[old_index] == 1);
                    grid_cell->border_zone = (fibrosis[old_index] == 2);
                    grid_cell->scar_type = tag[old_index];
                }
            } else {
                grid_cell->active = false;
            }
        }
        grid_cell = grid_cell->next;
    }

    fclose (file);

    // deallocate memory
    for (int l = 0; l < size; l++) {
        free (mesh_points[l]);
    }

    free (mesh_points);
    free (tag);
    free (fibrosis);
}

void set_cell_not_changeable (struct cell_node *c, double initialDiscretization) {

    double P1x, P1y, P1z;
    double P2x, P2y, P2z;
    double P3x, P3y, P3z;
    double P4x, P4y, P4z;
    double P5x, P5y, P5z;
    double P6x, P6y, P6z;
    double P7x, P7y, P7z;
    double P8x, P8y, P8z;
    double Cx, Cy, Cz;

    if (initialDiscretization == 100.0) {
        P1x = 6950;
        P1y = 50;
        P1z = 50;
        P2x = 6950;
        P2y = 19950;
        P2z = 50;
        P3x = 6950;
        P3y = 50;
        P3z = 2950;
        P4x = 6950;
        P4y = 19950;
        P4z = 2950;
        P5x = 50;
        P5y = 50;
        P5z = 50;
        P6x = 50;
        P6y = 19950;
        P6z = 50;
        P7x = 50;
        P7y = 50;
        P7z = 2950;
        P8x = 50;
        P8y = 19950;
        P8z = 2950;
        Cx = 3450;
        Cy = 9950;
        Cz = 1450;
    }

    else if (initialDiscretization == 200.0) {
        P1x = 6900;
        P1y = 100;
        P1z = 100;
        P2x = 6900;
        P2y = 19900;
        P2z = 100;
        P3x = 6900;
        P3y = 100;
        P3z = 2900;
        P4x = 6900;
        P4y = 19900;
        P4z = 2900;
        P5x = 100;
        P5y = 100;
        P5z = 100;
        P6x = 100;
        P6y = 19900;
        P6z = 100;
        P7x = 100;
        P7y = 100;
        P7z = 2900;
        P8x = 100;
        P8y = 19900;
        P8z = 2900;
        Cx = 3500;
        Cy = 9900;
        Cz = 1500;
    }

    else if (initialDiscretization == 125.0) {
        P1x = 6937.5;
        P1y = 62.5;
        P1z = 62.5;
        P2x = 6937.5;
        P2y = 19937.5;
        P2z = 62.5;
        P3x = 6937.5;
        P3y = 62.5;
        P3z = 2937.5;
        P4x = 6937.5;
        P4y = 19937.5;
        P4z = 2937.5;
        P5x = 62.5;
        P5y = 62.5;
        P5z = 62.5;
        P6x = 62.5;
        P6y = 19937.5;
        P6z = 62.5;
        P7x = 3937.5;
        P7y = 19937.5;
        P7z = 62.5;
        P8x = 62.5;
        P8y = 19937.5;
        P8z = 2937.5;
        Cx = 3437.5;
        Cy = 9937.5;
        Cz = 1562.5;
    }

    else if (initialDiscretization == 250.0) {
        P1x = 6875;
        P1y = 125;
        P1z = 125;
        P2x = 6875;
        P2y = 19875;
        P2z = 125;
        P3x = 6875;
        P3y = 125;
        P3z = 2875;
        P4x = 6875;
        P4y = 19875;
        P4z = 2875;
        P5x = 125;
        P5y = 125;
        P5z = 125;
        P6x = 125;
        P6y = 19875;
        P6z = 125;
        P7x = 125;
        P7y = 125;
        P7z = 2875;
        P8x = 125;
        P8y = 19875;
        P8z = 2875;
        Cx = 3375;
        Cy = 9875;
        Cz = 1125;
    }

    else {
        P1x = -1;
        P1y = -1;
        P1z = -1;
        P2x = -1;
        P2y = -1;
        P2z = -1;
        P3x = -1;
        P3y = -1;
        P3z = -1;
        P4x = -1;
        P4y = -1;
        P4z = -1;
        P5x = -1;
        P5y = -1;
        P5z = -1;
        P6x = -1;
        P6y = -1;
        P6z = -1;
        P7x = -1;
        P7y = -1;
        P7z = -1;
        P8x = -1;
        P8y = -1;
        P8z = -1;
        Cx = -1;
        Cy = -1;
        Cz = -1;
    }
    bool cannotChange = ((c->center_x == P1x) && (c->center_y == P1y) && (c->center_z == P1z));
    cannotChange |= ((c->center_x == P2x) && (c->center_y == P2y) && (c->center_z == P2z));
    cannotChange |= ((c->center_x == P3x) && (c->center_y == P3y) && (c->center_z == P3z));
    cannotChange |= ((c->center_x == P4x) && (c->center_y == P4y) && (c->center_z == P4z));
    cannotChange |= ((c->center_x == P5x) && (c->center_y == P5y) && (c->center_z == P5z));
    cannotChange |= ((c->center_x == P6x) && (c->center_y == P6y) && (c->center_z == P6z));
    cannotChange |= ((c->center_x == P7x) && (c->center_y == P7y) && (c->center_z == P7z));
    cannotChange |= ((c->center_x == P8x) && (c->center_y == P8y) && (c->center_z == P8z));
    cannotChange |= ((c->center_x == Cx) && (c->center_y == Cy) && (c->center_z == Cz));

    c->can_change = !cannotChange;
}

void set_plain_fibrosis (struct grid *the_grid, double phi, unsigned fib_seed) {

    print_to_stdout_and_file ("Making %.2lf %% of cells inactive\n", phi * 100.0);

    struct cell_node *grid_cell;

    if (fib_seed == 0)
        fib_seed = (unsigned)time (NULL) + getpid ();

    srand (fib_seed);

    print_to_stdout_and_file ("Using %u as seed\n", fib_seed);

    grid_cell = the_grid->first_cell;
    while (grid_cell != 0) {

        if (grid_cell->active) {
            double p = (double)(rand ()) / (RAND_MAX);
            if (p < phi)
                grid_cell->active = false;
            grid_cell->fibrotic = true;
        }
        grid_cell = grid_cell->next;
    }
}

void set_plain_sphere_fibrosis (struct grid *the_grid, double phi, double plain_center, double sphere_radius,
                                double bz_size, double bz_radius, unsigned fib_seed) {

    print_to_stdout_and_file ("Making %.2lf %% of cells inactive\n", phi * 100.0f);

    if (fib_seed == 0)
        fib_seed = (unsigned)time (NULL) + getpid ();

    srand (fib_seed);

    print_to_stdout_and_file ("Using %u as seed\n", fib_seed);

    double bz_radius_2 = pow (bz_radius, 2.0);
    double sphere_radius_2 = pow (sphere_radius, 2.0);

    struct cell_node *grid_cell;

    grid_cell = the_grid->first_cell;
    while (grid_cell != 0) {

        double distance = pow (grid_cell->center_x - plain_center, 2.0) + pow (grid_cell->center_y - plain_center, 2.0);

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
                double p = (double)(rand ()) / (RAND_MAX);
                if (p < phi)
                    grid_cell->active = false;
                grid_cell->can_change = false;
            } else if (grid_cell->border_zone) {
                double distance_from_center =
                        sqrt ((grid_cell->center_x - plain_center) * (grid_cell->center_x - plain_center) +
                              (grid_cell->center_y - plain_center) * (grid_cell->center_y - plain_center));
                distance_from_center = (distance_from_center - sphere_radius) / bz_size;
                double phi_local = phi - phi * distance_from_center;
                double p = (double)(rand ()) / (RAND_MAX);
                if (p < phi_local)
                    grid_cell->active = false;
                grid_cell->can_change = false;
            }
        }
        grid_cell = grid_cell->next;
    }
}

void set_human_mesh_fibrosis (struct grid *grid, double phi, unsigned seed, double big_scar_center_x,
                              double big_scar_center_y, double big_scar_center_z, double small_scar_center_x,
                              double small_scar_center_y, double small_scar_center_z) {

    if (seed == 0)
        seed = (unsigned) time(NULL) + getpid();

    srand(seed);

    print_to_stdout_and_file("Using %u as seed\n", seed);

    double bz_size_big = 0;
    double bz_size_small = 0;
    double dist_big = 0;
    double dist_small = 0;

    print_to_stdout_and_file("Calculating fibrosis using phi: %lf\n", phi);
    struct cell_node *gridCell = grid->first_cell;

    while (gridCell != NULL) {

        if (gridCell->active) {
            if (gridCell->fibrotic) {
                gridCell->can_change = false;
                double p = (double) (rand()) / (RAND_MAX);
                if (p < phi)
                    gridCell->active = false;
            } else if (gridCell->border_zone) {
                double centerX = gridCell->center_x;
                double centerY = gridCell->center_y;
                double centerZ = gridCell->center_z;
                if(gridCell->scar_type == 'b') {
                    dist_big = sqrt((centerX - big_scar_center_x) * (centerX - big_scar_center_x) +
                                    (centerY - big_scar_center_y) * (centerY - big_scar_center_y) +
                                    (centerZ - big_scar_center_z) * (centerZ - big_scar_center_z));
                    if (dist_big > bz_size_big) {
                        bz_size_big = dist_big;
                    }
                }
                else if(gridCell->scar_type == 's') {
                    dist_small = sqrt((centerX - small_scar_center_x) * (centerX - small_scar_center_x) +
                                      (centerY - small_scar_center_y) * (centerY - small_scar_center_y) +
                                      (centerZ - small_scar_center_z) * (centerZ - small_scar_center_z));
                    if (dist_small > bz_size_small) {
                        bz_size_small = dist_small;
                    }
                }
            }
        }
        gridCell = gridCell->next;
    }

    gridCell = grid->first_cell;
    while (gridCell != NULL) {

        if (gridCell->active) {
            if (gridCell->border_zone) {
                double centerX = gridCell->center_x;
                double centerY = gridCell->center_y;
                double centerZ = gridCell->center_z;
                if(gridCell->scar_type == 'b') {
                    dist_big = sqrt((centerX - big_scar_center_x) * (centerX - big_scar_center_x) +
                                    (centerY - big_scar_center_y) * (centerY - big_scar_center_y) +
                                    (centerZ - big_scar_center_z) * (centerZ - big_scar_center_z));
                    dist_big = dist_big / bz_size_big;
                    double phi_local = phi - phi * dist_big;

                    double p = (double) (rand()) / (RAND_MAX);
                    if (p < phi_local) {
                        gridCell->active = false;
                    }
                    gridCell->can_change = false;
                }
                else if(gridCell->scar_type == 's') {
                    dist_small = sqrt((centerX - small_scar_center_x) * (centerX - small_scar_center_x) +
                                      (centerY - small_scar_center_y) * (centerY - small_scar_center_y) +
                                      (centerZ - small_scar_center_z) * (centerZ - small_scar_center_z));
                    dist_small = dist_small / bz_size_small;
                    double phi_local = phi - phi * dist_small;

                    double p = (double) (rand()) / (RAND_MAX);
                    if (p < phi_local) {
                        gridCell->active = false;
                    }
                    gridCell->can_change = false;
                }
            }
        }
        gridCell = gridCell->next;
    }
}


void set_human_mesh_fibrosis_from_file(struct grid *grid, char type, const char *filename, int size) {


    FILE *file = fopen(filename, "r");

    if(!file) {
        printf("Error opening file %s!!\n", filename);
        exit(0);
    }


    double **scar_mesh = (double**) malloc(sizeof(double*)*size);
    for(int i=0; i< size; i++){
        scar_mesh[i] = (double*) malloc(sizeof(double) * 3);
        if(scar_mesh[i] == NULL) {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }
    double dummy1, dummy2; //unused values

    int i = 0;

    while ( !feof(file) ) {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf\n",&scar_mesh[i][0],&scar_mesh[i][1],&scar_mesh[i][2], &dummy1, &dummy2);
        i++;
    }

    fclose(file);

    sort_vector(scar_mesh, size);

    struct cell_node *grid_cell = grid->first_cell;
    while( grid_cell != 0 ) {

        double center_x = grid_cell->center_x;
        double center_y = grid_cell->center_y;
        double center_z = grid_cell->center_z;

        if( (grid_cell->face_length == 100.0) && (grid_cell->scar_type == type)) {
            int index = inside_mesh(scar_mesh, center_x, center_y, center_z, 0, size - 1);
            grid_cell->active = (index != -1);
        }

        grid_cell = grid_cell->next;
    }

    for(int k=0; k < size; k++){
        free(scar_mesh[k]);
    }

    free(scar_mesh);

}



