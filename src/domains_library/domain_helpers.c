//
// Created by sachetto on 19/10/17.
//

#include "domain_helpers.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"
#include "../utils/stop_watch.h"
#include "../utils/utils.h"
#include "mesh_info_data.h"

#include <float.h>
#include <math.h>
#include <time.h>

#include <unistd.h>

int set_cuboid_domain_mesh(struct grid *the_grid, real_cpu start_dx, real_cpu start_dy, real_cpu start_dz, real_cpu side_length_x, real_cpu side_length_y,
                           real_cpu side_length_z) {

    real_cpu real_side_length_x;
    real_cpu real_side_length_y;
    real_cpu real_side_length_z;

    int success = calculate_cuboid_side_lengths(start_dx, start_dy, start_dz, side_length_x, side_length_y, side_length_z, &real_side_length_x,
                                                &real_side_length_y, &real_side_length_z);

    if(!success) {
        return 0;
    }

    log_info("Initial mesh side length: %lf µm x %lf µm x %lf µm\n", real_side_length_x, real_side_length_y, real_side_length_z);
    log_info("Loading cuboid mesh with %lf µm x %lf µm x %lf µm using dx %lf µm, dy %lf µm, dz %lf µm\n", side_length_x, side_length_y, side_length_z, start_dx,
             start_dy, start_dz);

    int num_steps = get_num_refinement_steps_to_discretization(real_side_length_z, start_dz);

    initialize_and_construct_grid(the_grid, POINT3D(real_side_length_x, real_side_length_y, real_side_length_z));

    if((real_side_length_z / 2.0f) > side_length_z) {
        real_cpu aux = real_side_length_z / 2.0f;
        int remaining_refinements = num_steps;

        for(int i = 0; i < num_steps; i++) {
            set_cuboid_domain(the_grid, real_side_length_x, real_side_length_y, aux);
            refine_grid(the_grid, 1);
            derefine_grid_inactive_cells(the_grid);
            aux = aux / 2.0f;
            remaining_refinements--;
            if(aux <= side_length_z)
                break;
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

    the_grid->start_discretization.x = start_dx;
    the_grid->start_discretization.y = start_dy;
    the_grid->start_discretization.z = start_dz;

    return 1;
}

int set_square_mesh(struct config *config, struct grid *the_grid) {

    int num_layers = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, num_layers, config, "num_layers");

    real_cpu side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length, config, "side_length");

    real_cpu start_dx = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dx, config, "start_dx");

    real_cpu start_dy = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dy, config, "start_dy");

    real_cpu start_dz = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dz, config, "start_dz");

    return set_cuboid_domain_mesh(the_grid, start_dx, start_dy, start_dz, side_length, side_length, start_dz * num_layers);
}

uint32_t set_custom_mesh_from_file(struct grid *the_grid, const char *mesh_file, uint32_t num_volumes, double start_h, uint8_t num_extra_fields,
                                   set_custom_data_for_mesh_fn set_custom_data_for_mesh) {

    struct stop_watch sw = {0};
    start_stop_watch(&sw);

    FILE *file = fopen(mesh_file, "r");

    if(!file) {
        log_error_and_exit("Error opening mesh described in %s!!\n", mesh_file);
    }

    struct custom_mesh_basic_data_hash_entry *custom_mesh_data_hash = NULL;
    hmdefault(custom_mesh_data_hash, -1);

    real_cpu maxx = 0.0;
    real_cpu maxy = 0.0;
    real_cpu maxz = 0.0;
    real_cpu minx = DBL_MAX;
    real_cpu miny = DBL_MAX;
    real_cpu minz = DBL_MAX;

    bool load_custom_data = (set_custom_data_for_mesh != NULL && num_extra_fields > 0);

    real_cpu **custom_data = NULL;

    if(load_custom_data) {
        custom_data = MALLOC_ARRAY_OF_TYPE(real_cpu *, num_volumes);
        if(custom_data == NULL) {
            log_error_and_exit("Failed to allocate memory\n");
        }
        for(uint32_t i = 0; i < num_volumes; i++) {
            custom_data[i] = MALLOC_ARRAY_OF_TYPE(real_cpu, num_extra_fields);

            if(custom_data[i] == NULL) {
                log_error_and_exit("Failed to allocate memory\n");
            }
        }
    }

    char *line = NULL;
    size_t len;

    log_info("Start - reading mesh file\n");

    for(uint32_t i = 0; i < num_volumes; i++) {

        sds *data;
        int split_count;

        getline(&line, &len, file);

        char *tmp = line;
        data = sdssplit(tmp, ",", &split_count);

        if(split_count < 3) {
            log_error_and_exit("Not enough data to load the mesh geometry in line %d of file %s! [available=%d, required=3]\n", i + 1, mesh_file, split_count);
        }

        real_cpu cx = strtod(data[0], NULL);
        real_cpu cy = strtod(data[1], NULL);
        real_cpu cz = strtod(data[2], NULL);

        if(load_custom_data) {
            // indexes 3, 4 and 5 are not used in this function
            for(int d = 0; d < num_extra_fields; d++) {
                custom_data[i][d] = strtod(data[d + 6], NULL);
            }
        }

        hmput(custom_mesh_data_hash, POINT3D(cx, cy, cz), i);

        if(cx > maxx) {
            maxx = cx;
        }

        if(cx < minx) {
            minx = cx;
        }

        if(cy > maxy) {
            maxy = cy;
        }

        if(cy < miny) {
            miny = cy;
        }

        if(cz > maxz) {
            maxz = cz;
        }

        if(cz < minz) {
            minz = cz;
        }

        sdsfreesplitres(data, split_count);
    }

    log_info("Finish - reading mesh file\n");

    double cube_side = start_h;
    double min_cube_side = fmax(maxx, fmax(maxy, maxz)) + start_h;

    while(cube_side < min_cube_side) {
        cube_side = cube_side * 2;
    }

    double tmp_size = cube_side / 2.0;
    uint16_t num_ref = 0;

    while(tmp_size > start_h) {
        tmp_size = tmp_size / 2;
        num_ref++;
    }

    initialize_and_construct_grid(the_grid, SAME_POINT3D(cube_side));

    uint32_t num_loaded = 0;

    struct point_3d min_bound = POINT3D(minx - start_h, miny - start_h, minz - start_h);
    struct point_3d max_bound = POINT3D(maxx + start_h, maxy + start_h, maxz + start_h);

    log_info("\nStart - refining the initial cube (refining %d times with a cube side of %lf)\n", num_ref, cube_side);
    refine_grid_with_bounds(the_grid, num_ref, min_bound, max_bound);
    log_info("Finish - refining the initial cube\n\n");

    log_info("Loading grid with cube side of %lf\n", cube_side);

    FOR_EACH_CELL(the_grid) {
        real_cpu x = cell->center.x;
        real_cpu y = cell->center.y;
        real_cpu z = cell->center.z;

        if(x > maxx || y > maxy || z > maxz || x < minx || y < miny || z < minz) {
            cell->active = false;
        } else {

            struct point_3d p = POINT3D(x, y, z);
            int index = hmget(custom_mesh_data_hash, p);

            if(index != -1) {
                cell->active = true;
                cell->original_position_in_file = index;
                if(load_custom_data) {
                    set_custom_data_for_mesh(cell, custom_data[cell->original_position_in_file]);
                }

                num_loaded++;

            } else {
                cell->active = false;
            }
        }
    }

    if(num_loaded > 0) {
        the_grid->mesh_side_length.x = maxx + start_h;
        the_grid->mesh_side_length.y = maxy + start_h;
        the_grid->mesh_side_length.z = maxz + start_h;

        log_info("Cleaning grid\n");

        for(uint16_t r = 0; r < num_ref; r++) {
            derefine_grid_inactive_cells(the_grid);
        }
    }

    free(line);
    hmfree(custom_mesh_data_hash);

    if(custom_data) {
        for(uint32_t i = 0; i < num_volumes; i++) {
            free(custom_data[i]);
        }
        free(custom_data);
    }

    fclose(file);

    log_info("\nTime to load the mesh: %ld μs\n\n", stop_stop_watch(&sw));

    return num_loaded;
}

int calculate_cuboid_side_lengths(real_cpu start_dx, real_cpu start_dy, real_cpu start_dz, real_cpu side_length_x, real_cpu side_length_y,
                                  real_cpu side_length_z, real_cpu *real_side_length_x, real_cpu *real_side_length_y, real_cpu *real_side_length_z) {

    *real_side_length_x = start_dx * 2.0;
    *real_side_length_y = start_dy * 2.0;
    *real_side_length_z = start_dz * 2.0;

    real_cpu nx = side_length_x / start_dx;
    real_cpu ny = side_length_y / start_dy;
    real_cpu nz = side_length_z / start_dz;

    real_cpu proportion_dxdy = fmax(start_dx, start_dy) / fmin(start_dx, start_dy);
    real_cpu proportion_dxdz = fmax(start_dx, start_dz) / fmin(start_dx, start_dz);
    real_cpu proportion_dydz = fmax(start_dz, start_dy) / fmin(start_dz, start_dy);

    bool error = false;

    if(start_dx > start_dy) {
        if(side_length_x < side_length_y) {
            REPORT_ERROR_ON_FUNCTION_AND_CONTINUE("Incorrect configuration. If start_dx > start_dy, you need side_length_x > side_length_y");
            error = true;
        }
    }

    if(start_dx > start_dz) {
        if(side_length_x < side_length_z) {
            REPORT_ERROR_ON_FUNCTION_AND_CONTINUE("Incorrect configuration. If start_dx > start_dz, you need side_length_x > side_length_z");
            error = true;
        }
    }

    if(start_dy > start_dx) {
        if(side_length_y < side_length_x) {
            REPORT_ERROR_ON_FUNCTION_AND_CONTINUE("Incorrect configuration. If start_dy > start_dx, you need side_length_y > side_length_x");
            error = true;
        }
    }

    if(start_dy > start_dz) {
        if(side_length_y < side_length_z) {
            REPORT_ERROR_ON_FUNCTION_AND_CONTINUE("Incorrect configuration. If start_dy > start_dz, you need side_length_y > side_length_z");
            error = true;
        }
    }

    if(start_dz > start_dx) {
        if(side_length_z < side_length_x) {
            REPORT_ERROR_ON_FUNCTION_AND_CONTINUE("Incorrect configuration. If start_dz > start_dx, you need side_length_z > side_length_x");
            error = true;
        }
    }

    if(start_dz > start_dy) {
        if(side_length_z < side_length_y) {
            REPORT_ERROR_ON_FUNCTION_AND_CONTINUE("Incorrect configuration. If start_dz > start_dy, you need side_length_z > side_length_y");
            error = true;
        }
    }

    if(ceil(proportion_dxdy) != proportion_dxdy || ceil(proportion_dxdz) != proportion_dxdz || ceil(proportion_dydz) != proportion_dydz) {
        REPORT_ERROR_ON_FUNCTION_AND_CONTINUE("Incorrect configuration. start_dx, start_dy and start_dz need to be multiples");
        error = true;
    }

    if(ceil(nx) != nx) {
        sds error_str = sdscatprintf(sdsempty(), "start_dx: %lf is not multiple of side_length_x: %lf", start_dx, side_length_x);
        REPORT_ERROR_ON_FUNCTION_AND_CONTINUE(error_str);
        sdsfree(error_str);
        error = true;
    }
    if(ceil(ny) != ny) {
        sds error_str = sdscatprintf(sdsempty(), "start_dy: %lf is not multiple of side_length_y: %lf", start_dy, side_length_y);
        REPORT_ERROR_ON_FUNCTION_AND_CONTINUE(error_str);
        sdsfree(error_str);
        error = true;
    }
    if(ceil(nz) != nz) {
        sds error_str = sdscatprintf(sdsempty(), "start_dz: %lf is not multiple of side_length_z: %lf", start_dz, side_length_z);
        REPORT_ERROR_ON_FUNCTION_AND_CONTINUE(error_str);
        sdsfree(error_str);
        error = true;
    }

    if(error) {
        return 0;
    }

    while(*real_side_length_x < side_length_x) {
        *real_side_length_x *= 2.0;
    }

    while(*real_side_length_y < side_length_y) {
        *real_side_length_y *= 2.0;
    }

    while(*real_side_length_z < side_length_z) {
        *real_side_length_z *= 2.0;
    }

    int proportion_h;

    if(start_dx > start_dy) {
        proportion_h = (int)(start_dx / start_dy);
        if(*real_side_length_y >= *real_side_length_x || (*real_side_length_x / proportion_h) < *real_side_length_y) {
            *real_side_length_x = *real_side_length_y * proportion_h;
        } else {
            *real_side_length_y = *real_side_length_x / proportion_h;
            *real_side_length_z = *real_side_length_z / proportion_h;
        }

    } else if(start_dx < start_dy) {
        proportion_h = (int)(start_dy / start_dx);
        if(*real_side_length_x >= *real_side_length_y) {
            *real_side_length_y = *real_side_length_x * proportion_h;
            *real_side_length_z = *real_side_length_z * proportion_h;
        } else {
            *real_side_length_x = *real_side_length_y / proportion_h;
        }
    }

    if(start_dy > start_dz) {
        proportion_h = (int)(start_dy / start_dz);
        if(*real_side_length_z >= *real_side_length_y || *real_side_length_y / proportion_h < *real_side_length_z) {
            *real_side_length_y = *real_side_length_z * proportion_h;
            *real_side_length_x = *real_side_length_x * proportion_h;
        } else {
            *real_side_length_z = *real_side_length_y / proportion_h;
        }

    } else if(start_dy < start_dz) {
        proportion_h = (int)(start_dz / start_dy);
        if(*real_side_length_y > *real_side_length_z || *real_side_length_z / proportion_h < *real_side_length_y) {
            *real_side_length_z = *real_side_length_y * proportion_h;
        } else {
            *real_side_length_y = *real_side_length_z / proportion_h;
            *real_side_length_x = *real_side_length_x / proportion_h;
        }
    }

    if(start_dx == start_dy) {
        real_cpu aux = fmax(*real_side_length_x, *real_side_length_y);

        *real_side_length_x = aux;
        *real_side_length_y = aux;
    }

    if(start_dx == start_dz) {
        real_cpu aux = fmax(*real_side_length_x, *real_side_length_z);

        *real_side_length_x = aux;
        *real_side_length_z = aux;
    }

    if(start_dy == start_dz) {
        real_cpu aux = fmax(*real_side_length_y, *real_side_length_z);

        *real_side_length_y = aux;
        *real_side_length_z = aux;
    }

    return 1;
}

/**
 * Sets the current domain as a domain described in the N-version benchmark
 * (http://rsta.royalsocietypublishing.org/content/369/1954/4331)
 *
 */
void set_benchmark_domain(struct grid *the_grid) {
    struct cell_node *grid_cell = the_grid->first_cell;

    real_cpu sx, sy, sz;
    sx = 20000;
    sy = 7000;
    sz = 3000;

    while(grid_cell != 0) {
        grid_cell->active = (grid_cell->center.x < sx) && (grid_cell->center.y < sy) && (grid_cell->center.z < sz);
        grid_cell = grid_cell->next;
    }

    the_grid->mesh_side_length.x = sx;
    the_grid->mesh_side_length.y = sy;
    the_grid->mesh_side_length.z = sz;
}

void set_cuboid_domain(struct grid *the_grid, real_cpu size_x, real_cpu size_y, real_cpu size_z) {
    struct cell_node *grid_cell = the_grid->first_cell;

    while(grid_cell != 0) {
        grid_cell->active = (grid_cell->center.y < size_y) && (grid_cell->center.x < size_x) && (grid_cell->center.z < size_z);
        grid_cell = grid_cell->next;
    }

    the_grid->mesh_side_length.x = size_x;
    the_grid->mesh_side_length.y = size_y;
    the_grid->mesh_side_length.z = size_z;
}

void set_cell_not_changeable(struct cell_node *c, real_cpu initialDiscretization) {

    real_cpu P1x, P1y, P1z;
    real_cpu P2x, P2y, P2z;
    real_cpu P3x, P3y, P3z;
    real_cpu P4x, P4y, P4z;
    real_cpu P5x, P5y, P5z;
    real_cpu P6x, P6y, P6z;
    real_cpu P7x, P7y, P7z;
    real_cpu P8x, P8y, P8z;
    real_cpu Cx, Cy, Cz;

    if(initialDiscretization == 100.0) {
        P1x = 50;
        P1y = 6950;
        P1z = 50;

        P2x = 19950;
        P2y = 6950;
        P2z = 50;

        P3x = 50;
        P3y = 6950;
        P3z = 2950;

        P4x = 19950;
        P4y = 6950;
        P4z = 2950;

        P5x = 50;
        P5y = 50;
        P5z = 50;

        P6x = 19950;
        P6y = 50;
        P6z = 50;

        P7x = 50;
        P7y = 50;
        P7z = 2950;

        P8x = 19950;
        P8y = 50;
        P8z = 2950;

        Cx = 9950;
        Cy = 3450;
        Cz = 1450;
    }

    else if(initialDiscretization == 200.0) {
        P1x = 100;
        P1y = 6900;
        P1z = 100;

        P2x= 19900;
        P2y = 6900;
        P2z = 100;

        P3x = 100;
        P3y = 6900;
        P3z = 2900;

        P4x = 19900;
        P4y = 6900;
        P4z = 2900;

        P5x = 100;
        P5y = 100;
        P5z = 100;

        P6x = 19900;
        P6y = 100;
        P6z = 100;

        P7x = 100;
        P7y = 100;
        P7z = 2900;

        P8x = 19900;
        P8y = 100;
        P8z = 2900;

        Cx = 9900;
        Cy = 3500;
        Cz = 1500;
    }

    else if(initialDiscretization == 125.0) {
        P1x = 62.5;
        P1y = 6937.5;
        P1z = 62.5;

        P2x = 19937.5;
        P2y = 6937.5;
        P2z = 62.5;

        P3x = 62.5;
        P3y = 6937.5;
        P3z = 2937.5;

        P4x = 19937.5;
        P4y = 6937.5;
        P4z = 2937.5;

        P5x = 62.5;
        P5y = 62.5;
        P5z = 62.5;

        P6x = 19937.5;
        P6y = 62.5;
        P6z = 62.5;

        P7x = 19937.5;
        P7y = 3937.5;
        P7z = 62.5;

        P8x = 19937.5;
        P8y = 62.5;
        P8z = 2937.5;

        Cx = 9937.5;
        Cy = 3437.5;
        Cz = 1562.5;
    }

    else if(initialDiscretization == 250.0) {
        P1x = 125;
        P1y = 6875;
        P1z = 125;

        P2x = 19875;
        P2y = 6875;
        P2z = 125;

        P3x = 125;
        P3y = 6875;
        P3z = 2875;

        P4x = 19875;
        P4y = 6875;
        P4z = 2875;

        P5x = 125;
        P5y = 125;
        P5z = 125;

        P6x = 19875;
        P6y = 125;
        P6z = 125;

        P7x = 125;
        P7y = 125;
        P7z = 2875;

        P8x = 19875;
        P8y = 125;
        P8z = 2875;

        Cx = 9875;
        Cy = 3375;
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
    bool cannotChange = ((c->center.x == P1x) && (c->center.y == P1y) && (c->center.z == P1z));
    cannotChange |= ((c->center.x == P2x) && (c->center.y == P2y) && (c->center.z == P2z));
    cannotChange |= ((c->center.x == P3x) && (c->center.y == P3y) && (c->center.z == P3z));
    cannotChange |= ((c->center.x == P4x) && (c->center.y == P4y) && (c->center.z == P4z));
    cannotChange |= ((c->center.x == P5x) && (c->center.y == P5y) && (c->center.z == P5z));
    cannotChange |= ((c->center.x == P6x) && (c->center.y == P6y) && (c->center.z == P6z));
    cannotChange |= ((c->center.x == P7x) && (c->center.y == P7y) && (c->center.z == P7z));
    cannotChange |= ((c->center.x == P8x) && (c->center.y == P8y) && (c->center.z == P8z));
    cannotChange |= ((c->center.x == Cx) && (c->center.y == Cy) && (c->center.z == Cz));

    c->can_change = !cannotChange;

    if(cannotChange) {
        printf("Cannot change %lf, %lf, %lf\n", c->center.x, c->center.y, c->center.z);
    }
}

void set_plain_fibrosis(struct grid *the_grid, real_cpu phi, unsigned fib_seed) {

    log_info("Making %.2lf %% of cells inactive\n", phi * 100.0);

    log_info("Using %u as seed\n", fib_seed);

    struct cell_node *grid_cell;

    grid_cell = the_grid->first_cell;
    while(grid_cell != 0) {

        if(grid_cell->active) {
            real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
            if(p < phi) {
                grid_cell->active = false;
            }

            INITIALIZE_FIBROTIC_INFO(grid_cell);
            FIBROTIC(grid_cell) = true;
        }
        grid_cell = grid_cell->next;
    }
}

void set_plain_source_sink_fibrosis(struct grid *the_grid, real_cpu channel_width, real_cpu channel_length) {

    log_info("Making upper and down left corner inactive !\n");

    bool inside;

    //   real_cpu side_length_x = the_grid->mesh_side_length.x;
    real_cpu side_length_y = the_grid->mesh_side_length.y;
    //    real_cpu side_length_z = the_grid->mesh_side_length.z;

    real_cpu region_height = (side_length_y - channel_width) / 2.0;

    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;

    while(grid_cell != 0) {

        if(grid_cell->active) {

            real_cpu x = grid_cell->center.x;
            real_cpu y = grid_cell->center.y;
            //            real_cpu z = grid_cell->center.z;

            // Check region 1
            inside = (x >= 0.0) && (x <= channel_length) && (y >= 0.0) && (y <= region_height);

            // Check region 2
            inside |= (x >= 0.0) && (x <= channel_length) && (y >= region_height + channel_width) && (y <= side_length_y);

            if(inside) {
                grid_cell->active = false;
            }

            INITIALIZE_FIBROTIC_INFO(grid_cell);
            FIBROTIC(grid_cell) = true;
        }
        grid_cell = grid_cell->next;
    }
}

void set_plain_sphere_fibrosis(struct grid *the_grid, real_cpu phi, real_cpu plain_center, real_cpu sphere_radius, real_cpu bz_size, real_cpu bz_radius,
                               unsigned fib_seed) {

    log_info("Making %.2lf %% of cells inactive\n", phi * 100.0f);

    if(fib_seed == 0)
        fib_seed = (unsigned)time(NULL) + getpid();

    srand(fib_seed);

    log_info("Using %u as seed\n", fib_seed);

    real_cpu bz_radius_2 = pow(bz_radius, 2.0);
    real_cpu sphere_radius_2 = pow(sphere_radius, 2.0);

    struct cell_node *grid_cell;

    grid_cell = the_grid->first_cell;
    while(grid_cell != 0) {

        real_cpu distance = pow(grid_cell->center.x - plain_center, 2.0) + pow(grid_cell->center.y - plain_center, 2.0);

        if(grid_cell->active) {

            INITIALIZE_FIBROTIC_INFO(grid_cell);

            if(distance <= bz_radius_2) {
                if(distance <= sphere_radius_2) {
                    FIBROTIC(grid_cell) = true;
                } else {
                    BORDER_ZONE(grid_cell) = true;
                }
            }
        }
        grid_cell = grid_cell->next;
    }

    grid_cell = the_grid->first_cell;

    while(grid_cell != 0) {

        if(grid_cell->active) {
            if(FIBROTIC(grid_cell)) {
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                if(p < phi)
                    grid_cell->active = false;
                grid_cell->can_change = false;
            } else if(BORDER_ZONE(grid_cell)) {
                real_cpu distance_from_center = sqrt((grid_cell->center.x - plain_center) * (grid_cell->center.x - plain_center) +
                                                     (grid_cell->center.y - plain_center) * (grid_cell->center.y - plain_center));
                distance_from_center = (distance_from_center - sphere_radius) / bz_size;
                real_cpu phi_local = phi - phi * distance_from_center;
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                if(p < phi_local)
                    grid_cell->active = false;
                grid_cell->can_change = false;
            }
        }
        grid_cell = grid_cell->next;
    }
}

void set_cube_sphere_fibrosis(struct grid *the_grid, real_cpu phi, real_cpu sphere_center[3], real_cpu sphere_radius, unsigned fib_seed) {

    log_info("Making %.2lf %% of cells inactive\n", phi * 100.0f);

    if(fib_seed == 0)
        fib_seed = (unsigned)time(NULL) + getpid();

    srand(fib_seed);

    log_info("Using %u as seed\n", fib_seed);

    real_cpu sphere_radius_2 = pow(sphere_radius, 2.0);

    struct cell_node *grid_cell;

    grid_cell = the_grid->first_cell;
    while(grid_cell != 0) {

        real_cpu distance = pow(grid_cell->center.x - sphere_center[0], 2.0) + pow(grid_cell->center.y - sphere_center[1], 2.0) + pow(grid_cell->center.z - sphere_center[2], 2.0);

        if(grid_cell->active) {

            INITIALIZE_FIBROTIC_INFO(grid_cell);
            if(distance <= sphere_radius_2) {
                FIBROTIC(grid_cell) = true;
            }
        }
        grid_cell = grid_cell->next;
    }

    grid_cell = the_grid->first_cell;

    while(grid_cell != 0) {

        if(grid_cell->active) {
            if(FIBROTIC(grid_cell)) {
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                if(p < phi)
                    grid_cell->active = false;
                grid_cell->can_change = false;
            } 
        }
        grid_cell = grid_cell->next;
    }
}

void set_plain_sphere_fibrosis_without_inactivating(struct grid *the_grid, real_cpu plain_center, real_cpu sphere_radius, real_cpu bz_radius) {

    real_cpu bz_radius_2 = pow(bz_radius, 2.0);
    real_cpu sphere_radius_2 = pow(sphere_radius, 2.0);

    struct cell_node *grid_cell;

    grid_cell = the_grid->first_cell;
    while(grid_cell != 0) {

        real_cpu distance = pow(grid_cell->center.x - plain_center, 2.0) + pow(grid_cell->center.y - plain_center, 2.0);

        if(grid_cell->active) {

            INITIALIZE_FIBROTIC_INFO(grid_cell);

            if(distance <= bz_radius_2) {
                if(distance <= sphere_radius_2) {
                    FIBROTIC(grid_cell) = true;
                } else {
                    BORDER_ZONE(grid_cell) = true;
                }
            }
        }
        grid_cell = grid_cell->next;
    }
}

void set_fibrosis_from_file(struct grid *grid, const char *filename, int size) {

    FILE *file = fopen(filename, "r");

    if(!file) {
        printf("Error opening file %s!!\n", filename);
        exit(0);
    }

    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * size);

    for(int i = 0; i < size; i++) {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 7);
        if(scar_mesh[i] == NULL) {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    for(int i = 0; i < size; i++) {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3], &scar_mesh[i][4],
               &scar_mesh[i][5], &scar_mesh[i][6]);
    }

    fclose(file);

    OMP(parallel for)
    for(int j = 0; j < size; j++) {

        struct cell_node *grid_cell = grid->first_cell;

        real_cpu b_center_x = scar_mesh[j][0];
        real_cpu b_center_y = scar_mesh[j][1];

        real_cpu b_h_dx = scar_mesh[j][3];
        real_cpu b_h_dy = scar_mesh[j][4];

        bool active = (bool)(scar_mesh[j][6]);

        int c = 0;
        while(grid_cell != 0) {

            if(grid_cell->active) {

                real_cpu center_x = grid_cell->center.x;
                real_cpu center_y = grid_cell->center.y;

                real_cpu half_dy = grid_cell->discretization.y / 2.0;

                if(FIBROTIC_INFO(grid_cell) == NULL) {
                    INITIALIZE_FIBROTIC_INFO(grid_cell);
                    FIBROTIC(grid_cell) = 1;
                }

                struct point_3d p;

                p.x = b_center_x + b_h_dx;
                p.y = b_center_y - b_h_dy;

                if(center_x == b_center_x && center_y + half_dy <= p.x && center_y - half_dy >= p.y) {
                    grid_cell->active = active;
                    c++;
                }
            }
            if(c == 4)
                break;

            grid_cell = grid_cell->next;
        }
    }

    for(int k = 0; k < size; k++) {
        free(scar_mesh[k]);
    }

    free(scar_mesh);
}

void set_plain_fibrosis_inside_region(struct grid *the_grid, real_cpu phi, unsigned fib_seed, const double min_x, const double max_x, const double min_y,
                                      const double max_y, const double min_z, const double max_z) {

    log_info("Making %.2lf %% of cells inside the region inactive\n", phi * 100.0);

    struct cell_node *grid_cell;

    if(fib_seed == 0)
        fib_seed = (unsigned)time(NULL) + getpid();

    srand(fib_seed);

    log_info("Using %u as seed\n", fib_seed);

    grid_cell = the_grid->first_cell;
    while(grid_cell != 0) {
        real center_x = grid_cell->center.x;
        real center_y = grid_cell->center.y;
        real center_z = grid_cell->center.z;

        if(center_x >= min_x && center_x <= max_x && center_y >= min_y && center_y <= max_y && center_z >= min_z && center_z <= max_z) {
            if(grid_cell->active) {
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                if(p < phi) {
                    grid_cell->active = false;
                }

                INITIALIZE_FIBROTIC_INFO(grid_cell);
                FIBROTIC(grid_cell) = true;
            }
        }

        grid_cell = grid_cell->next;
    }
}

int calc_num_refs(real_cpu start_h, real_cpu desired_h) {
    int num_refs = 0;
    while(start_h > desired_h) {
        start_h = start_h/2.0;
        num_refs++;
    }

    return num_refs;
}
