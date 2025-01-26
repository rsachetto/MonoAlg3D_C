//
// Created by sachetto on 01/10/17.
//

#include <unistd.h>

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../domains_library/mesh_info_data.h"
#include "helper_functions.h"

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

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_fibrosis));

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

    real_cpu dx = the_grid->start_discretization.x;
    real_cpu dy = the_grid->start_discretization.y;
    real_cpu dz = the_grid->start_discretization.z;

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

// 'libten_tusscher_tt3_mixed_endo_mid_epi.so' with transmurality and fibrosis (all cells healthy)
SET_EXTRA_DATA (set_extra_data_mixed_tt3) {
    uint32_t num_active_cells = the_grid->num_active_cells;
    real side_length = the_grid->mesh_side_length.x;
    struct cell_node ** ac = the_grid->active_cells;

    //
    struct extra_data_for_tt3 *extra_data = NULL;
    extra_data = set_common_tt3_data(config, num_active_cells);

    // Divide the domain in three sections (ENDO/MID/EPI)
    // The percentages were taken from the ToRORd paper (Transmural experiment)
    real side_length_endo = side_length*0.45;
    real side_length_mid = side_length_endo + side_length*0.25;
    real side_length_epi = side_length_mid + side_length*0.3;
	
	int i;

    // Transmurality and fibrosis tags
	OMP(parallel for)
    for (int i = 0; i < num_active_cells; i++) {

        real center_x = ac[i]->center.x;

        // Tag the model transmurality
        // ENDO=0, MID=1, EPI=2
        if (center_x < side_length_endo)
            extra_data->transmurality[i] = 0.0;
        else if (center_x >= side_length_endo && center_x < side_length_mid)
            extra_data->transmurality[i] = 1.0;
        else
            extra_data->transmurality[i] = 2.0;

        // Tag the fibrosis region
        extra_data->fibrosis[i] = 1.0;
    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_tt3));

    return (void*)extra_data;
}

// Initial condition - 'libToRORd_fkatp_mixed_endo_mid_epi.so' + transmurality + current modifiers (plain and cuboid)
SET_EXTRA_DATA(set_extra_data_mixed_torord_fkatp_epi_mid_endo) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    real side_length = the_grid->mesh_side_length.x;
    struct cell_node ** ac = the_grid->active_cells;

    // The percentages were taken from the ToRORd paper (Transmural experiment)
    real side_length_endo = side_length*0.45;
    real side_length_mid = side_length_endo + side_length*0.25;
    real side_length_epi = side_length_mid + side_length*0.3;

    struct extra_data_for_torord *extra_data = NULL;
    extra_data = set_common_torord_data(config, num_active_cells);

    OMP(parallel for)
    for (int i = 0; i < num_active_cells; i++) {

        real center_x = ac[i]->center.x;

        // ENDO
        if (center_x < side_length_endo)
            extra_data->transmurality[i] = 0.0;
        // MID
        else if (center_x >= side_length_endo && center_x < side_length_mid)
            extra_data->transmurality[i] = 1.0;
        // EPI
        else
            extra_data->transmurality[i] = 2.0;
    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_torord));

    return (void*)extra_data;

}

// Initial condition - 'libToRORd_dynCl_mixed_endo_mid_epi.so' + transmurality + current modifiers (plain and cuboid)
SET_EXTRA_DATA(set_extra_data_mixed_torord_dynCl_epi_mid_endo) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    real side_length = the_grid->mesh_side_length.x;
    struct cell_node ** ac = the_grid->active_cells;

    // The percentages were taken from the ToRORd paper (Transmural experiment)
    real side_length_endo = side_length*0.45;
    real side_length_mid = side_length_endo + side_length*0.25;
    real side_length_epi = side_length_mid + side_length*0.3;

    struct extra_data_for_torord *extra_data = NULL;
    extra_data = set_common_torord_dyncl_data(config, num_active_cells);

    OMP(parallel for)
    for (int i = 0; i < num_active_cells; i++) {

        real center_x = ac[i]->center.x;

        // ENDO
        if (center_x < side_length_endo)
            extra_data->transmurality[i] = 0.0;
        // MID
        else if (center_x >= side_length_endo && center_x < side_length_mid)
            extra_data->transmurality[i] = 1.0;
        // EPI
        else
            extra_data->transmurality[i] = 2.0;
    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_torord));

    return (void*)extra_data;
}

// Initial condition - 'libToRORd_Land_mixed_endo_mid_epi.so' + transmurality + current modifiers (plain and cuboid)
SET_EXTRA_DATA(set_extra_data_mixed_torord_Land_epi_mid_endo) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    real side_length = the_grid->mesh_side_length.x;
    struct cell_node ** ac = the_grid->active_cells;

    // The percentages were taken from the ToRORd paper (Transmural experiment)
    real side_length_endo = side_length*0.45;
    real side_length_mid = side_length_endo + side_length*0.25;
    real side_length_epi = side_length_mid + side_length*0.3;

    struct extra_data_for_torord_land *extra_data = NULL;
    extra_data = set_common_torord_Land_data(config, num_active_cells);

    OMP(parallel for)
    for (int i = 0; i < num_active_cells; i++) {

        real center_x = ac[i]->center.x;

        // ENDO
        if (center_x < side_length_endo)
            extra_data->transmurality[i] = 0.0;
        // MID
        else if (center_x >= side_length_endo && center_x < side_length_mid)
            extra_data->transmurality[i] = 1.0;
        // EPI
        else
            extra_data->transmurality[i] = 2.0;
    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_torord_land));

    return (void*)extra_data;

}

// Initial condition - 'libToRORd_Land_mixed_endo_mid_epi.so' + current modifiers (cell, plain and cuboid)
SET_EXTRA_DATA(set_extra_data_mixed_torord_Land_same_celltype) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node ** ac = the_grid->active_cells;

    struct extra_data_for_torord_land *extra_data = NULL;
    extra_data = set_common_torord_Land_data(config, num_active_cells);

    // All cells will be the same type
    OMP(parallel for)
    for (int i = 0; i < num_active_cells; i++) {
        // ENDO
        extra_data->transmurality[i] = 0.0;
        // MID
        //extra_data->transmurality[i] = 1.0;
        // EPI
        //extra_data->transmurality[i] = 2.0;
    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_torord_land));

    return (void*)extra_data;

}

// Current modifiers - 'libtrovato_2020.so'
SET_EXTRA_DATA(set_extra_data_trovato) {

    if (!the_grid->purkinje) {
        fprintf(stderr,"[ERR] There is no Purkinje network configured for this mesh!\n");
        exit(EXIT_FAILURE);
    }

    uint32_t num_purkinje_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **purkinje_ac = the_grid->purkinje->purkinje_cells;

    struct extra_data_for_trovato *extra_data = NULL;

    extra_data = set_common_trovato_data(config, num_purkinje_active_cells);

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_trovato));

    return (void*)extra_data;

}
