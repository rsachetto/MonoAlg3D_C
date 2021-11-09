//
// Created by sachetto on 01/10/17.
//

#include <unistd.h>

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../domains_library/mesh_info_data.h"
#include "helper_functions.h"


SET_EXTRA_DATA(set_extra_data_for_fibrosis_sphere) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node ** ac = the_grid->active_cells;

    struct extra_data_for_fibrosis *extra_data = NULL;

    real plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_center, config, "plain_center");

    real border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, border_zone_size, config, "border_zone_size");

    real sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sphere_radius, config, "sphere_radius");

    extra_data = set_common_schemia_data(config, num_active_cells);


    OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++) {

        if(FIBROTIC(ac[i])) {
            extra_data->fibrosis[i] = 0.0;
        }
        else if(BORDER_ZONE(ac[i])) {

            real center_x = (real)ac[i]->center.x;
            real center_y = (real)ac[i]->center.y;
            //TODO: Maybe we want the distance from the Z as well
            //real center_z = (real)ac[i]->center_z;

            real distanceFromCenter = sqrtf((center_x - plain_center)*(center_x - plain_center) + (center_y - plain_center)*(center_y - plain_center));
            distanceFromCenter = (distanceFromCenter - sphere_radius)/border_zone_size;
            extra_data->fibrosis[i] = distanceFromCenter;

        }
        else {
            extra_data->fibrosis[i] = 1.0;
        }

    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_fibrosis));

    return (void*)extra_data;
}

SET_EXTRA_DATA(set_extra_data_for_fibrosis_plain) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    struct extra_data_for_fibrosis *extra_data = NULL;

    extra_data = set_common_schemia_data(config, num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        extra_data->fibrosis[i] = 0.0;
    }

   SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_fibrosis));

    return (void*)extra_data;
}

SET_EXTRA_DATA(set_extra_data_for_no_fibrosis) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    struct extra_data_for_fibrosis *extra_data = NULL;

    extra_data = set_common_schemia_data(config, num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        extra_data->fibrosis[i] = 1.0;
    }

    return (void*)extra_data;
}

SET_EXTRA_DATA(set_extra_data_for_benchmark) {

    *extra_data_size = sizeof(real)*19;

    real *initial_conditions = (real*)malloc(*extra_data_size);

    // Initial conditions  // Var      Units          Initial value
    initial_conditions[ 0] = -85.423f;  // V;       millivolt;     -85.423
    initial_conditions[ 1] = 0.0165;   // Xr1;     dimensionless; 0.0165
    initial_conditions[ 2] = 0.473;    // Xr2;     dimensionless; 0.473
    initial_conditions[ 3] = 0.0174;   // Xs;      dimensionless; 0.0174
    initial_conditions[ 4] = 0.00165;  // m;       dimensionless; 0.00165
    initial_conditions[ 5] = 0.749;    // h;       dimensionless; 0.749
    initial_conditions[ 6] = 0.6788;   // j;       dimensionless; 0.6788
    initial_conditions[ 7] = 3.288e-5; // d;       dimensionless; 3.288e-5
    initial_conditions[ 8] = 0.7026;   // f;       dimensionless; 0.7026
    initial_conditions[ 9] = 0.9526;   // f2;      dimensionless; 0.9526
    initial_conditions[10] = 0.9942;   // fCass;   dimensionless; 0.9942
    initial_conditions[11] = 0.999998; // s;       dimensionless; 0.999998
    initial_conditions[12] = 2.347e-8; // r;       dimensionless; 2.347e-8
    initial_conditions[13] = 0.000153; // Ca_i;    millimolar;    0.000153
    initial_conditions[14] = 4.272;    // Ca_SR;   millimolar;    4.272
    initial_conditions[15] = 0.00042;  // Ca_ss;   millimolar;    0.00042
    initial_conditions[16] = 0.8978;   // R_prime; dimensionless; 0.8978
    initial_conditions[17] = 10.132;   // Na_i;    millimolar;    10.132
    initial_conditions[18] = 138.52;   // K_i;     millimolar;    138.52

    return (void*)initial_conditions;
}

SET_EXTRA_DATA (set_mixed_model_if_x_less_than)
{
    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells);

    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    real x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, x_limit, config, "x_limit");

    int i;
    bool inside;

    OMP(parallel for)
    for (i = 0; i < num_active_cells; i++)
    {
        real center_x = ac[i]->center.x;

        inside = (center_x <= x_limit);

        if (inside)
            mapping[i] = 0;
        else
            mapping[i] = 1;
    }

    return (void*)mapping;
}

SET_EXTRA_DATA (set_mixed_model_purkinje_and_tissue)
{
    uint32_t num_active_tissue_cells = the_grid->num_active_cells;
    uint32_t num_active_purkinje_cells = the_grid->purkinje->num_active_purkinje_cells;
    uint32_t num_active_cells = num_active_tissue_cells + num_active_purkinje_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells + 1);

    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    int i;

    // Purkinje section
    OMP(parallel for)
    for (i = 0; i < num_active_purkinje_cells; i++) {
        mapping[i] = 0;
    }

    // Tissue section
    OMP(parallel for)
    for (i = num_active_purkinje_cells; i < num_active_cells; i++) {
        mapping[i] = 1;
    }

    return (void*)mapping;
}

SET_EXTRA_DATA(set_extra_data_mixed_model_epi_mid_endo) {

    uint32_t num_eq = 43;   // ToRORd number of equations
    uint32_t num_active_cells = the_grid->num_active_cells;
    real side_length = the_grid->mesh_side_length.x;

    // The percentages were taken from the ToRORd paper (Transmural experiment)
    real side_length_endo = side_length*0.45;
    real side_length_mid = side_length_endo + side_length*0.25;
    real side_length_epi = side_length_mid + side_length*0.3;

    // The extra data size is the initial state vector for each celltype plus the mapping of each cell
    *extra_data_size = sizeof(real)*(num_active_cells) + sizeof(real)*num_eq*3;

    real *extra_data = (real*)malloc(*extra_data_size);

    // Set the initial conditions (celltype=ENDO)
    int offset = 0;
    extra_data[0] = -88.7638;
    extra_data[1] = 0.0111;
    extra_data[2] = 7.0305e-5;
    extra_data[3] = 12.1025;
    extra_data[4] = 12.1029;
    extra_data[5] = 142.3002;
    extra_data[6] = 142.3002;
    extra_data[7] = 1.5211;
    extra_data[8] = 1.5214;
    extra_data[9] = 8.1583e-05;
    extra_data[10] = 8.0572e-4;
    extra_data[11] = 0.8286;
    extra_data[12] = 0.8284;
    extra_data[13] = 0.6707;
    extra_data[14] = 0.8281;
    extra_data[15] = 1.629e-4;
    extra_data[16] = 0.5255;
    extra_data[17] = 0.2872;
    extra_data[18] = 9.5098e-4;
    extra_data[19] = 0.9996;
    extra_data[20] = 0.5936;
    extra_data[21] = 4.8454e-4;
    extra_data[22] = 0.9996;
    extra_data[23] = 0.6538;
    extra_data[24] = 8.1084e-9;
    extra_data[25] = 1.0;
    extra_data[26] = 0.939;
    extra_data[27] = 1.0;
    extra_data[28] = 0.9999;
    extra_data[29] = 1.0;
    extra_data[30] = 1.0;
    extra_data[31] = 1.0;
    extra_data[32] = 6.6462e-4;
    extra_data[33] = 0.0012;
    extra_data[34] = 7.0344e-4;
    extra_data[35] = 8.5109e-4;
    extra_data[36] = 0.9981;
    extra_data[37] = 1.3289e-5;
    extra_data[38] = 3.7585e-4;
    extra_data[39] = 0.248;
    extra_data[40] = 1.7707e-4;
    extra_data[41] = 1.6129e-22;
    extra_data[42] = 1.2475e-20;

    // Set the initial conditions (celltype=EPI)
    offset = num_eq;
    extra_data[0+offset] = -89.1400;
    extra_data[1+offset] = 0.0129;
    extra_data[2+offset] = 5.7672e-05;
    extra_data[3+offset] = 12.8363;
    extra_data[4+offset] = 12.8366;
    extra_data[5+offset] = 142.6951;
    extra_data[6+offset] = 142.6951;
    extra_data[7+offset] = 1.8119;
    extra_data[8+offset] = 1.8102;
    extra_data[9+offset] = 6.6309e-05;
    extra_data[10+offset] = 0.00074303;
    extra_data[11+offset] = 0.8360;
    extra_data[12+offset] = 0.8359;
    extra_data[13+offset] = 0.6828;
    extra_data[14+offset] = 0.8357;
    extra_data[15+offset] = 0.00015166;
    extra_data[16+offset] = 0.5401;
    extra_data[17+offset] = 0.3034;
    extra_data[18+offset] = 0.00092716;
    extra_data[19+offset] = 0.9996;
    extra_data[20+offset] = 0.9996;
    extra_data[21+offset] = 0.0004724;
    extra_data[22+offset] = 0.9996;
    extra_data[23+offset] = 0.9996;
    extra_data[24+offset] = 0;
    extra_data[25+offset] = 1.0;
    extra_data[26+offset] = 0.9485;
    extra_data[27+offset] = 1.0;
    extra_data[28+offset] = 0.9999;
    extra_data[29+offset] = 1.0;
    extra_data[30+offset] = 1.0;
    extra_data[31+offset] = 1.0;
    extra_data[32+offset] = 0.00030853;
    extra_data[33+offset] = 0.00053006;
    extra_data[34+offset] = 0.00067941;
    extra_data[35+offset] = 0.00082869;
    extra_data[36+offset] = 0.9982;
    extra_data[37+offset] = 9.5416e-06;
    extra_data[38+offset] = 0.00027561;
    extra_data[39+offset] = 0.2309;
    extra_data[40+offset] = 0.00016975;
    extra_data[41+offset] = 2.8189e-24;
    extra_data[42+offset] = 0;

    // Set the initial conditions (celltype=MID)
    offset = num_eq*2;
    extra_data[0+offset] = -89.1704;
    extra_data[1+offset] = 0.0192;
    extra_data[2+offset] = 6.5781e-05;
    extra_data[3+offset] = 15.0038;
    extra_data[4+offset] = 15.0043;
    extra_data[5+offset] = 143.0403;
    extra_data[6+offset] = 143.0402;
    extra_data[7+offset] = 1.9557;
    extra_data[8+offset] = 1.9593;
    extra_data[9+offset] = 8.166e-05;
    extra_data[10+offset] = 0.00073818;
    extra_data[11+offset] = 0.8365;
    extra_data[12+offset] = 0.8363;
    extra_data[13+offset] = 0.6838;
    extra_data[14+offset] = 0.8358;
    extra_data[15+offset] = 0.00015079;
    extra_data[16+offset] = 0.5327;
    extra_data[17+offset] = 0.2834;
    extra_data[18+offset] = 0.00092527;
    extra_data[19+offset] = 0.9996;
    extra_data[20+offset] = 0.5671;
    extra_data[21+offset] = 0.00047143;
    extra_data[22+offset] = 0.9996;
    extra_data[23+offset] = 0.6261;
    extra_data[24+offset] = 0;
    extra_data[25+offset] = 1.0;
    extra_data[26+offset] = 0.92;
    extra_data[27+offset] = 1.0;
    extra_data[28+offset] = 0.9998;
    extra_data[29+offset] = 1.0;
    extra_data[30+offset] = 1.0;
    extra_data[31+offset] = 1.0;
    extra_data[32+offset] = 0.00051399;
    extra_data[33+offset] = 0.0012;
    extra_data[34+offset] = 0.00069560;
    extra_data[35+offset] = 0.00082672;
    extra_data[36+offset] = 0.9979;
    extra_data[37+offset] = 1.8784e-05;
    extra_data[38+offset] = 0.00054206;
    extra_data[39+offset] = 0.2653;
    extra_data[40+offset] = 0.00016921;
    extra_data[41+offset] = 0;
    extra_data[42+offset] = 0;

    offset = num_eq*3;

    struct cell_node ** ac = the_grid->active_cells;

    int i;

    OMP(parallel for)
    for (i = 0; i < num_active_cells; i++) {

        real center_x = ac[i]->center.x;

        // ENDO
        if (center_x < side_length_endo)
            extra_data[i+offset] = 0.0;
        // MID
        else if (center_x >= side_length_endo && center_x < side_length_mid)
            extra_data[i+offset] = 2.0;
        // EPI
        else
            extra_data[i+offset] = 1.0;

    }

    return (void*)extra_data;

}

struct uv {
    real_cpu u;
    real_cpu v;
};

struct original_volumes_hash_entry {
    struct point_3d key;
    struct uv value;
};

struct cells_hash_entry {
    struct point_3d key;
    struct cell_node *value;
};


struct point_3d find_next_point(struct original_volumes_hash_entry *original_volumes, int l, struct point_3d p, struct point_3d *n, int np) {

    double dist;
    double min_dist = 100000.0;
    int min_idx = -1;

    int continue_loop = false;
    struct point_3d aux;

    for(int i = 0; i < l; i++) {

       aux = original_volumes[i].key;
       double ax = aux.x;
       double ay = aux.y;

        for(int j = 0; j < np; j++) {
            if(ax == n[j].x && ay == n[j].y) {
                continue_loop = true;
                break;
            }
        }

        if(continue_loop) {
            continue_loop = false;
            continue;
        }

        double a = ax - p.x;
        double b = ay - p.y;

        dist = a*a + b*b;

        if(dist < min_dist) {
            min_dist = dist;
            min_idx = i;
        }
    }

    if(min_idx > -1) {
        return original_volumes[min_idx].key;
    }
    else {
        return POINT3D(0,0,0);
    }

}

SET_EXTRA_DATA(set_extra_data_for_spiral_fhn) {


    size_t num_sv_entries = 2;
    uint32_t num_active_cells = the_grid->num_active_cells;


    *extra_data_size = num_active_cells * num_sv_entries;

    real *sv_cpu;
    sv_cpu = CALLOC_ARRAY_OF_TYPE(real, *extra_data_size);

    char *file_u;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_u, config, "file_u");

    char *file_v;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_v, config, "file_v");

    bool interpolate = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(interpolate, config, "interpolate");

    FILE *file_U = fopen(file_u, "r");
    FILE *file_V = fopen(file_v, "r");

    double value;

    struct cell_node *cell = the_grid->first_cell;

    real_cpu dx = 100;
    real_cpu dy = 100;
    real_cpu dz = 100;

    int nx = 18000/dx;
    int ny = 18000/dy;

    real_cpu half_dx = dx/2.0;
    real_cpu half_dy = dy/2.0;
    real_cpu half_dz = dz/2.0;

    real_cpu center_x = 0;
    real_cpu center_y = 50;
    real_cpu center_z = 50;


    struct cells_hash_entry *volumes = NULL;

    FOR_EACH_CELL(the_grid) {
        if(cell->active)
            hmput(volumes, cell->center, cell);
    }

    if(interpolate) {

        printf("-----Interpolate-----\n");

        struct original_volumes_hash_entry *original_volumes = NULL;

        for(int y = 0; y < 180; y++) {

            center_x = 50;

            for(int x = 0; x < 180; x++) {

                struct point_3d center = POINT3D(center_x, center_y, center_z);
                struct uv og;

                fscanf(file_U,"%le",&value);
                og.u = value;
                fscanf(file_V,"%le",&value);
                og.v = value;

                hmput(original_volumes, center, og);

                center_x += 100;
            }
            center_y += 100;
        }

        center_y = half_dy;
        center_z = half_dz;

        int l = 180*180;

        for(int y = 0; y < ny; y++) {
            center_x = half_dx;

            for(int x = 0; x < nx; x++) {

                struct point_3d *not_include = NULL;

                struct point_3d p1  = find_next_point(original_volumes, l, POINT3D(center_x, center_y, 50), not_include, 0);
                arrput(not_include, p1);

                struct point_3d p2  = find_next_point(original_volumes, l, POINT3D(center_x, center_y, 50), not_include, 1);
                arrput(not_include, p2);

                struct point_3d p3  = find_next_point(original_volumes, l, POINT3D(center_x, center_y, 50), not_include, 2);
                arrput(not_include, p3);

                struct point_3d p4  = find_next_point(original_volumes, l, POINT3D(center_x, center_y, 50), not_include, 3);
                arrfree(not_include);

                int idx_b = hmgeti(original_volumes, p1);
                int idx_f = hmgeti(original_volumes, p2);
                int idx_u = hmgeti(original_volumes, p3);
                int idx_d = hmgeti(original_volumes, p4);

                real_cpu u_b = 0;
                real_cpu u_f = 0;
                real_cpu u_u = 0;
                real_cpu u_d = 0;

                real_cpu v_b = 0;
                real_cpu v_f = 0;
                real_cpu v_u = 0;
                real_cpu v_d = 0;


                //printf("%lf, %lf, %lf\n", center_x, center_y, center_z);
                //printf("%lf, %lf, %lf, %d\n", p1.x, p1.y, p1.z, idx_b);
                //printf("%lf, %lf, %lf, %d\n", p2.x, p2.y, p2.z, idx_f);

                //printf("%lf, %lf, %lf, %d\n", p3.x, p3.y, p3.z, idx_u);
                //printf("%lf, %lf, %lf, %d\n", p4.x, p4.y, p4.z, idx_d);


                if(idx_b != -1) {
                    u_b = original_volumes[idx_b].value.u;
                    v_b = original_volumes[idx_b].value.v;
                }

                if(idx_f != -1) {
                    u_f = original_volumes[idx_f].value.u;
                    v_f = original_volumes[idx_f].value.v;
                }


                if(idx_u != -1) {
                    u_u = original_volumes[idx_u].value.u;
                    v_u = original_volumes[idx_u].value.v;
                }

                if(idx_d != -1) {
                    u_d = original_volumes[idx_d].value.u;
                    v_d = original_volumes[idx_d].value.v;
                }

                real_cpu u = (u_b + u_f + u_u + u_d)/4.0;
                real_cpu v = (v_b + v_f + v_u + v_d)/4.0;


                struct point_3d center = POINT3D(center_x, center_y, center_z);
                struct cell_node *c = hmget(volumes, center);

                sv_cpu[c->sv_position] = u;

                sv_cpu[num_active_cells + c->sv_position] = v;

                center_x += dx;
            }
            center_y += dy;
        }

    } else {
        for(int y = 0; y < ny; y++) {

            center_x = half_dx;

            for(int x = 0; x < nx; x++) {

                struct point_3d center = POINT3D(center_x, center_y, center_z);
                struct cell_node *c = hmget(volumes, center);

                fscanf(file_U,"%le",&value);
                sv_cpu[c->sv_position] = value;

                fscanf(file_V,"%le",&value);
                sv_cpu[num_active_cells + c->sv_position] = value;

                center_x += dx;
            }

            center_y += dy;
        }
    }
    return sv_cpu;
}
