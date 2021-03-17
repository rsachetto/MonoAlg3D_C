#include <ctype.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "3dparty/sds/sds.h"
#include "3dparty/stb_ds.h"
#include "common_types/common_types.h"
#include "config/config_parser.h"
#include "utils/file_utils.h"

struct elem {
    struct point_3d n1;
    struct point_3d n2;
    struct point_3d n3;
    struct point_3d n4;
};

int closest_node(struct elem *elements, struct point_3d alg_node) {

    int n = arrlen(elements);
    int nt = omp_get_max_threads();

    double min_local[nt];
    int index_local[nt];

    for(int i = 0; i < nt; i++) {
        min_local[i] = FLT_MAX;
    }

    int index = 0;

    OMP(parallel for)
    for(int element_index = 0; element_index < n; element_index++) {

        int tid = omp_get_thread_num();

        double dist;
        double x[4] = {elements[element_index].n1.x, elements[element_index].n2.x, elements[element_index].n3.x, elements[element_index].n4.x};
        double y[4] = {elements[element_index].n1.y, elements[element_index].n2.y, elements[element_index].n3.y, elements[element_index].n4.y};
        double z[4] = {elements[element_index].n1.z, elements[element_index].n2.z, elements[element_index].n3.z, elements[element_index].n4.z};

        struct point_3d centroid;

        centroid.x = (x[0] + x[1] + x[2] + x[3]) / 4.0;
        centroid.y = (y[0] + y[1] + y[2] + y[3]) / 4.0;
        centroid.z = (z[0] + z[1] + z[2] + z[3]) / 4.0;

        dist = sqrt(pow(alg_node.x - centroid.x, 2) + pow(alg_node.y - centroid.y, 2) + pow(alg_node.z - centroid.z, 2));

        if(dist < min_local[tid]) {
            min_local[tid] = dist;
            index_local[tid] = element_index;
        }
    }

    double min = FLT_MAX;
    for(int i = 0; i < nt; i++) {
        if(min_local[i] < min) {
            min = min_local[i];
            index = index_local[i];
        }
    }

    return index;
}

void convert_fibers(struct fibers_conversion_options *options) {

    FILE *fibers_file = open_file_or_exit(options->fibers_file, "r");
    FILE *elements_file = open_file_or_exit(options->ele_file, "r");
    FILE *nodes_file = open_file_or_exit(options->nodes_file, "r");
    FILE *alg_file = open_file_or_exit(options->alg_file, "r");
    FILE *outfile = open_file_or_exit(options->output_file, "w");

    point3d_array nodes = NULL;
    struct elem *elements = NULL;
    double *dists = NULL;
    struct fiber_coords *alg_fiber_coords = NULL;
    struct fiber_coords *fiber_coords = NULL;

    char *line = NULL;
    size_t len;

    sds *points;
    int split_count;

    double x;
    double y;
    double z;

    double max_x_new = 77.25;
    double max_x_old = 143.525;
    double min_x_old = 65.8178;

    double max_y_new = 89.25;
    double min_y_old = 28.4905;
    double max_y_old = 117.842;

    double max_z_new = 99.25;
    double min_z_old = 366.488;
    double max_z_old = 465.862;

    printf("Reading nodes\n");
    getline(&line, &len, nodes_file); // skip first line
    while((getline(&line, &len, nodes_file)) != -1) {
        char *tmp = line;

        while(!isblank(*tmp))
            tmp++;
        while(isblank(*tmp))
            tmp++;

        points = sdssplit(tmp, " ", &split_count);

        x = strtod(points[0], NULL);
        y = strtod(points[1], NULL);
        z = strtod(points[2], NULL);

        double new_x = (((x - min_x_old) * max_x_new) / (max_x_old - min_x_old)) * 1000.0;
        double new_y = (((y - min_y_old) * max_y_new) / (max_y_old - min_y_old)) * 1000.0;
        double new_z = (((z - min_z_old) * max_z_new) / (max_z_old - min_z_old)) * 1000.0;

        arrput(nodes, POINT3D(new_x, new_y, new_z));
        sdsfreesplitres(points, split_count);
    }
    free(line);
    line = NULL;

    printf("Reading elements\n");
    getline(&line, &len, elements_file); // skip first line
    while((getline(&line, &len, elements_file)) != -1) {

        char *tmp = line;
        while(!isblank(*tmp))
            tmp++;
        while(isblank(*tmp))
            tmp++;

        points = sdssplit(tmp, "\t", &split_count);

        int n1 = strtol(points[0], NULL, 10) - 1;
        int n2 = strtol(points[1], NULL, 10) - 1;
        int n3 = strtol(points[2], NULL, 10) - 1;
        int n4 = strtol(points[3], NULL, 10) - 1;

        struct elem e;
        e.n1 = nodes[n1];
        e.n2 = nodes[n2];
        e.n3 = nodes[n3];
        e.n4 = nodes[n4];
        arrput(elements, e);
        arrput(dists, -1.0);

        sdsfreesplitres(points, split_count);
    }
    free(line);
    line = NULL;

    arrsetcap(fiber_coords, arrlen(elements));

    printf("Reading fibers\n");
    getline(&line, &len, fibers_file); // skip first line

    while((getline(&line, &len, fibers_file)) != -1) {
        char *tmp = line;

        while(!isblank(*tmp))
            tmp++;
        while(isblank(*tmp))
            tmp++;

        points = sdssplit(tmp, " ", &split_count);
        struct fiber_coords f_coords;

        f_coords.f[0] = strtod(points[0], NULL);
        f_coords.f[1] = strtod(points[1], NULL);
        f_coords.f[2] = strtod(points[2], NULL);

        f_coords.s[0] = strtod(points[3], NULL);
        f_coords.s[1] = strtod(points[4], NULL);
        f_coords.s[2] = strtod(points[5], NULL);

        f_coords.n[0] = strtod(points[6], NULL);
        f_coords.n[1] = strtod(points[7], NULL);
        f_coords.n[2] = strtod(points[8], NULL);

        arrput(fiber_coords, f_coords);
        sdsfreesplitres(points, split_count);
    }
    free(line);
    line = NULL;

    printf("Reading ALG mesh\n");

    int count = 0;
    while((getline(&line, &len, alg_file)) != -1) {
        points = sdssplit(line, ",", &split_count);
        x = strtod(points[0], NULL);
        y = strtod(points[1], NULL);
        z = strtod(points[2], NULL);

        sdsfreesplitres(points, split_count);
        struct point_3d alg_center = POINT3D(x, y, z);
        int closest_index = closest_node(elements, alg_center);

        arrput(alg_fiber_coords, fiber_coords[closest_index]);
        count++;
        if(count % 2500 == 0) {
            printf("%d lines processed\n", count);
        }
    }

    free(line);

    for(int i = 0; i < arrlen(alg_fiber_coords); i++) {
        fprintf(outfile, "%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", alg_fiber_coords[i].f[0], alg_fiber_coords[i].f[1],
                alg_fiber_coords[i].f[2], alg_fiber_coords[i].s[0], alg_fiber_coords[i].s[1], alg_fiber_coords[i].s[2], alg_fiber_coords[i].n[0],
                alg_fiber_coords[i].n[1], alg_fiber_coords[i].n[2]);
    }

    fclose(outfile);
}
int main(int argc, char **argv) {

    struct fibers_conversion_options *options = new_fibers_conversion_options();

    parse_fibers_conversion_options(argc, argv, options);

    if(!(options->fibers_file && options->ele_file && options->nodes_file && options->alg_file)) {
        display_fibers_conversion_usage(argv);
        return EXIT_FAILURE;
    }

    printf("Fibers file: %s\n", options->fibers_file);
    printf("Elements file: %s\n", options->ele_file);
    printf("Nodes file file: %s\n", options->nodes_file);
    printf("Alg mesh file file: %s\n", options->alg_file);
    printf("Output file: %s\n\n", options->output_file);

    convert_fibers(options);

    return EXIT_SUCCESS;
}
