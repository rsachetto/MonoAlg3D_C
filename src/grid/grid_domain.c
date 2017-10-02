//
// Created by sachetto on 01/10/17.
//

#include "../utils/point_hash.h"
#include "../utils/utils.h"
#include "grid.h"
#include <float.h>
#include <math.h>

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

void set_plain_domain (struct grid *the_grid, float sizeX, float sizeY, float sizeZ) {
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

    // TODO: point_hash map for C
    // TPoint3DMap myMap;
    struct point_hash *p_hash = hash_create ();
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

        ////TODO: hash map for C
        point3d.x = a[i][0];
        point3d.y = a[i][1];
        point3d.z = a[i][2];

        // myMap[TPoint3D(a[i])] = fibrosis;
        hash_insert (p_hash, point3d, fibrosis);

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

        // TODO: point_hash map for C
        // grid_cell->fibrotic = (myMap[TPoint3D(x,y,z)] == 1);
        // grid_cell->border_zone = (myMap[TPoint3D(x,y,z)] == 2);
        grid_cell->fibrotic = (hash_search (p_hash, point3d) == 1);
        grid_cell->border_zone = (hash_search (p_hash, point3d) == 2);

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
    hash_destroy (p_hash);
}

// TODO: da pra mudar essa funcao toda usando o hashmap
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

    // TODO: @Incomplete: point_hash map for C
    // TPoint3DMap myMap;
    struct point_hash *p_hash = hash_create ();
    struct point_3d point3d;

    while (i < size) {
        fscanf (file, "%lf,%lf,%lf,%lf,%d\n", &a[i][0], &a[i][1], &a[i][2], &b, &fibrosis);
        // TODO: @Incomplete: point_hash map for C
        //  myMap[TPoint3D(a[i])] = fibrosis;
        point3d.x = a[i][0];
        point3d.y = a[i][1];
        point3d.z = a[i][2];

        // myMap[TPoint3D(a[i])] = fibrosis;
        hash_insert (p_hash, point3d, fibrosis);

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

        // TODO: @Incomplete: point_hash map for C
        // grid_cell->fibrotic = (myMap[TPoint3D(x,y,z)] == 1);
        // grid_cell->border_zone = (myMap[TPoint3D(x,y,z)] == 2);
        grid_cell->fibrotic = (hash_search (p_hash, point3d) == 1);
        grid_cell->border_zone = (hash_search (p_hash, point3d) == 2);

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
    hash_destroy (p_hash);
}

void set_human_mesh (struct grid *the_grid, const char *file_name) {
    set_custom_mesh (the_grid, file_name, 2025252);
}

void set_human_sub_mesh (struct grid *the_grid, const char *file_name, double minx, double maxx,
                         double miny, double maxy, double minz, double maxz) {
    set_custom_mesh_with_bounds (the_grid, file_name, 2025252, minx, maxx, miny, maxy, minz, maxz);
}

void set_rabbit_mesh (struct grid *the_grid, const char *mesh_file) {
    set_custom_mesh (the_grid, mesh_file, 470197);
}

void initialize_grid_with_mouse_mesh (struct grid *the_grid, const char *mesh_file) {

    if (the_grid == NULL) {
        fprintf(stderr, "set_mouse_mesh: the_grid is NULL. Exiting!");
        exit(10);
    }

    initialize_and_construct_grid(the_grid, 6400.0);

    refine_grid(the_grid, 5);

    printf("Loading Mouse Heart Mesh\n");;
    set_custom_mesh(the_grid, mesh_file, 96195);

    printf("Cleaning grid\n");

    int i;
    for( i = 0; i < 5; i++ ) {
         derefine_grid_inactive_cells(the_grid);
    }

}

void save_grid_domain (struct grid *the_grid, const char *file_name) {
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
