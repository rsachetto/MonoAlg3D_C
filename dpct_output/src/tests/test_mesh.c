////
//// Created by sachetto on 06/10/17.
////
#include <criterion/criterion.h>

#include "../3dparty/ini_parser/ini.h"
#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../alg/grid/grid.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"

int test_cuboid_mesh(char *start_dx, char *start_dy, char *start_dz, char *side_length_x, char *side_length_y,
                     char *side_length_z, bool save, bool compress, bool binary, int id) {

    struct grid *grid = new_grid();
    struct config *domain_config;

    domain_config = alloc_and_init_config_data();

    shput_dup_value(domain_config->config_data, "start_dx", start_dx);
    shput_dup_value(domain_config->config_data, "start_dy", start_dy);
    shput_dup_value(domain_config->config_data, "start_dz", start_dz);

    domain_config->main_function_name = strdup("initialize_grid_with_cuboid_mesh");
    shput_dup_value(domain_config->config_data, "name", "Test cuboid");

    shput(domain_config->config_data, "side_length_x", strdup(side_length_x));
    shput(domain_config->config_data, "side_length_y", strdup(side_length_y));
    shput(domain_config->config_data, "side_length_z", strdup(side_length_z));

    init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

    int success = ((set_spatial_domain_fn *)domain_config->main_function)(domain_config, grid);

    if(!success) {
        clean_and_free_grid(grid);
        free_config_data(domain_config);
        return 0;
    }

    order_grid_cells(grid);

    real_cpu sx = grid->mesh_side_length.x;
    real_cpu sy = grid->mesh_side_length.y;
    real_cpu sz = grid->mesh_side_length.z;

    real_cpu start_dx_value = strtod(start_dx, NULL);
    real_cpu start_dy_value = strtod(start_dy, NULL);
    real_cpu start_dz_value = strtod(start_dz, NULL);

    real_cpu nx = sx / start_dx_value;
    real_cpu ny = sy / start_dy_value;
    real_cpu nz = sz / start_dz_value;

    struct cell_node *cell = grid->first_cell;

    real_cpu max_x = 0;
    real_cpu max_y = 0;
    real_cpu max_z = 0;

    while(cell) {

        if(cell->active) {
            if(cell->center.x > max_x) {
                max_x = cell->center.x;
            }

            if(cell->center.y > max_y) {
                max_y = cell->center.y;
            }

            if(cell->center.z > max_z) {
                max_z = cell->center.z;
            }
        }

        cell = cell->next;
    }

    if(save) {
        struct config *save_mesh_config = alloc_and_init_config_data();

        save_mesh_config->init_function_name = strdup("init_save_as_vtk_or_vtu");
        save_mesh_config->main_function_name = strdup("save_as_vtu");
        save_mesh_config->end_function_name = strdup("end_save_as_vtk_or_vtu");
        shput_dup_value(save_mesh_config->config_data, "output_dir", "./tests_bin");
        shput_dup_value(save_mesh_config->config_data, "print_rate", "1");

        sds file_prefix = sdscatprintf(sdsempty(), "test_%s_%s_%s_%s_%s_%s_%d", start_dx, start_dy, start_dz,
                                       side_length_x, side_length_y, side_length_z, id);

        init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_result");

        shput(save_mesh_config->config_data, "file_prefix", strdup(file_prefix));

        if(compress) {
            shput(save_mesh_config->config_data, "compress", strdup("yes"));
        } else if(binary) {
            shput(save_mesh_config->config_data, "binary", strdup("yes"));
        }

        shput(save_mesh_config->config_data, "save_pvd", strdup("no"));

        struct time_info ti = ZERO_TIME_INFO;

        ((init_save_mesh_fn *)save_mesh_config->init_function)(save_mesh_config);
        ((save_mesh_fn *)save_mesh_config->main_function)(&ti, save_mesh_config, grid, NULL, NULL);
        ((end_save_mesh_fn *)save_mesh_config->end_function)(save_mesh_config, grid);

        free_config_data(save_mesh_config);
        sdsfree(file_prefix);
    }

    cr_assert_float_eq(max_x + (start_dx_value / 2.0), atof(side_length_x), 1e-16);
    cr_assert_float_eq(max_y + (start_dy_value / 2.0), atof(side_length_y), 1e-16);
    cr_assert_float_eq(max_z + (start_dz_value / 2.0), atof(side_length_z), 1e-16);
    cr_assert_eq(nx * ny * nz, grid->num_active_cells);

    clean_and_free_grid(grid);
    free_config_data(domain_config);

    return 1;
}

int compare_two_binary_files(FILE *fp1, FILE *fp2) {
    char ch1, ch2;
    int flag = 0;

    while(((ch1 = fgetc(fp1)) != EOF) && ((ch2 = fgetc(fp2)) != EOF)) {
        /*
         * character by character comparision
         * if equal then continue by comparing till the end of files
         */
        if(ch1 == ch2) {
            flag = 1;
            continue;
        }
        /*
         * If not equal then returns the byte position
         */
        else {
            fseek(fp1, -1, SEEK_CUR);
            flag = 0;
            break;
        }
    }

    if(flag == 0) {
        return ftell(fp1) + 1;
    } else {
        return -1;
    }
}

Test(mesh_load, cuboid_mesh_100_100_100_1000_1000_1000) {
    int success = test_cuboid_mesh("100", "100", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test(mesh_load, cuboid_mesh_200_100_100_1000_1000_1000) {
    int success = test_cuboid_mesh("200", "100", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test(mesh_load, cuboid_mesh_100_200_100_1000_1000_1000) {

    int success = test_cuboid_mesh("100", "200", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test(mesh_load, cuboid_mesh_100_100_200_1000_1000_1000) {
    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test(mesh_load, cuboid_mesh_100_100_100_1000_1000_2000) {
    int success = test_cuboid_mesh("100", "100", "100", "1000", "1000", "2000", false, false, false, 0);
    cr_assert(success);
}

Test(mesh_load, cuboid_mesh_150_150_150_1500_1500_1500) {
    int success = test_cuboid_mesh("150", "150", "150", "1500", "1500", "1500", false, false, false, 0);
    cr_assert(success);
}

Test(mesh_load, cuboid_mesh_150_150_150_1500_1500_3000) {
    int success = test_cuboid_mesh("150", "150", "150", "1500", "1500", "3000", false, false, false, 0);
    cr_assert(success);
}
Test(mesh_load, cuboid_mesh_300_150_150_1500_1500_3000) {
    int success = test_cuboid_mesh("300", "150", "150", "1500", "1500", "3000", false, false, false, 0);
    cr_assert(!success);
}

Test(mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_binary) {

    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", true, false, true, 1);
    cr_assert(success);

    FILE *f1 = fopen("tests_bin/test_100_100_200_1000_1000_1000_1_it_0.vtu", "r");
    FILE *f2 = fopen("tests_bin/gold_vtu_mesh_binary.vtu", "r");

    success = compare_two_binary_files(f1, f2);

    fclose(f1);
    fclose(f2);

    cr_assert(success == -1);
}

Test(mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_compressed) {

    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", true, true, false, 2);
    cr_assert(success);

    FILE *f1 = fopen("tests_bin/test_100_100_200_1000_1000_1000_2_it_0.vtu", "r");
    FILE *f2 = fopen("tests_bin/gold_vtu_mesh_compressed.vtu", "r");

    cr_assert(f1);
    cr_assert(f2);

    success = compare_two_binary_files(f1, f2);

    fclose(f1);
    fclose(f2);

    cr_assert(success == -1);
}

Test(mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_plain) {

    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", true, false, false, 3);
    cr_assert(success);

    FILE *f1 = fopen("tests_bin/test_100_100_200_1000_1000_1000_3_it_0.vtu", "r");
    FILE *f2 = fopen("tests_bin/gold_vtu_mesh.vtu", "r");

    cr_assert(f1);
    cr_assert(f2);

    success = compare_two_binary_files(f1, f2);

    fclose(f1);
    fclose(f2);

    cr_assert(success == -1);
}

#define assert_node(direction, expected_neighbour)                                                                     \
    node = get_cell_neighbour_with_same_refinement_level(cell, direction);                                             \
    cr_assert_eq(node, expected_neighbour);

#define assert_node_center(direction, expect_active)                                                                   \
    node = get_cell_neighbour_with_same_refinement_level(cell, direction);                                             \
    if(node) {                                                                                                         \
        cr_assert_eq(node->active, (expect_active));                                                                   \
        cr_assert(CENTER_EQUALS(node->center, (direction##_CELL)));                                                    \
    } else {                                                                                                           \
        cr_assert(!(expect_active));                                                                                   \
    }

#define FRONT_CELL TRANSLATE(center, 0, 0, d)
#define BACK_CELL TRANSLATE(center, 0, 0, -d)

#define TOP_CELL TRANSLATE(center, 0, d, 0)
#define DOWN_CELL TRANSLATE(center, 0, -d, 0)

#define RIGHT_CELL TRANSLATE(center, d, 0, 0)
#define LEFT_CELL TRANSLATE(center, -d, 0, 0)

#define FRONT_TOP_CELL TRANSLATE(center, 0, d, d)
#define FRONT_DOWN_CELL TRANSLATE(center, 0, -d, d)
#define FRONT_RIGHT_CELL TRANSLATE(center, d, 0, d)
#define FRONT_LEFT_CELL TRANSLATE(center, -d, 0, d)
#define FRONT_TOP_RIGHT_CELL TRANSLATE(center, d, d, d)
#define FRONT_TOP_LEFT_CELL TRANSLATE(center, -d, d, d)
#define FRONT_DOWN_RIGHT_CELL TRANSLATE(center, d, -d, d)
#define FRONT_DOWN_LEFT_CELL TRANSLATE(center, -d, -d, d)
#define BACK_TOP_CELL TRANSLATE(center, 0, d, -d)
#define BACK_DOWN_CELL TRANSLATE(center, 0, -d, -d)
#define BACK_RIGHT_CELL TRANSLATE(center, d, 0, -d)
#define BACK_LEFT_CELL TRANSLATE(center, -d, 0, -d)
#define BACK_TOP_RIGHT_CELL TRANSLATE(center, d, d, -d)
#define BACK_TOP_LEFT_CELL TRANSLATE(center, -d, d, -d)
#define BACK_DOWN_RIGHT_CELL TRANSLATE(center, d, -d, -d)
#define BACK_DOWN_LEFT_CELL TRANSLATE(center, -d, -d, -d)
#define TOP_RIGHT_CELL TRANSLATE(center, d, d, 0)
#define TOP_LEFT_CELL TRANSLATE(center, -d, d, 0)
#define DOWN_RIGHT_CELL TRANSLATE(center, d, -d, 0)
#define DOWN_LEFT_CELL TRANSLATE(center, -d, -d, 0)

Test(cell_conectors, start_cube) {
    struct grid *grid = new_grid();

    initialize_and_construct_grid(grid, POINT3D(600, 600, 600));
    order_grid_cells(grid);

    // right_front_top_cell is the first cell
    struct cell_node *right_front_top_cell = grid->first_cell;
    struct cell_node *right_front_down_cell = right_front_top_cell->neighbours[DOWN];
    struct cell_node *right_back_top_cell = right_front_top_cell->neighbours[BACK];
    struct cell_node *right_back_down_cell = right_back_top_cell->neighbours[DOWN];
    struct cell_node *left_front_top_cell = right_front_top_cell->neighbours[LEFT];
    struct cell_node *left_front_down_cell = right_front_down_cell->neighbours[LEFT];
    struct cell_node *left_back_top_cell = right_back_top_cell->neighbours[LEFT];
    struct cell_node *left_back_down_cell = right_back_down_cell->neighbours[LEFT];

    FOR_EACH_CELL(grid) {

        struct cell_node *node = NULL;

        if(CENTER_EQUALS(cell->center, POINT3D(150, 150, 150))) {

            assert_node(FRONT, left_front_down_cell);
            assert_node(BACK, NULL);
            assert_node(TOP, left_back_top_cell);
            assert_node(DOWN, NULL);
            assert_node(RIGHT, right_back_down_cell);
            assert_node(LEFT, NULL);
            assert_node(FRONT_TOP, left_front_top_cell);
            assert_node(FRONT_DOWN, NULL);
            assert_node(FRONT_RIGHT, right_front_down_cell);
            assert_node(FRONT_LEFT, NULL);
            assert_node(FRONT_TOP_RIGHT, right_front_top_cell);
            assert_node(FRONT_TOP_LEFT, NULL);
            assert_node(FRONT_DOWN_RIGHT, NULL);
            assert_node(FRONT_DOWN_LEFT, NULL);
            assert_node(BACK_TOP, NULL);
            assert_node(BACK_DOWN, NULL);
            assert_node(BACK_RIGHT, NULL);
            assert_node(BACK_LEFT, NULL);
            assert_node(BACK_TOP_RIGHT, NULL);
            assert_node(BACK_TOP_LEFT, NULL);
            assert_node(BACK_DOWN_RIGHT, NULL);
            assert_node(BACK_DOWN_LEFT, NULL);
            assert_node(TOP_RIGHT, right_back_top_cell);
            assert_node(TOP_LEFT, NULL);
            assert_node(DOWN_RIGHT, NULL);
            assert_node(DOWN_LEFT, NULL);
        }

        if(CENTER_EQUALS(cell->center, POINT3D(150, 150, 450))) {

            assert_node(FRONT, NULL);
            assert_node(BACK, left_back_down_cell);
            assert_node(TOP, left_front_top_cell);
            assert_node(DOWN, NULL);
            assert_node(RIGHT, right_front_down_cell);
            assert_node(LEFT, NULL);
            assert_node(FRONT_TOP, NULL);
            assert_node(FRONT_DOWN, NULL);
            assert_node(FRONT_RIGHT, NULL);
            assert_node(FRONT_LEFT, NULL);
            assert_node(FRONT_TOP_RIGHT, NULL);
            assert_node(FRONT_TOP_LEFT, NULL);
            assert_node(FRONT_DOWN_RIGHT, NULL);
            assert_node(FRONT_DOWN_LEFT, NULL);
            assert_node(BACK_TOP, left_back_top_cell);
            assert_node(BACK_DOWN, NULL);
            assert_node(BACK_RIGHT, right_back_down_cell);
            assert_node(BACK_LEFT, NULL);
            assert_node(BACK_TOP_RIGHT, right_back_top_cell);
            assert_node(BACK_TOP_LEFT, NULL);
            assert_node(BACK_DOWN_RIGHT, NULL);
            assert_node(BACK_DOWN_LEFT, NULL);
            assert_node(TOP_RIGHT, right_front_top_cell);
            assert_node(TOP_LEFT, NULL);
            assert_node(DOWN_RIGHT, NULL);
            assert_node(DOWN_LEFT, NULL);
        }

        if(CENTER_EQUALS(cell->center, POINT3D(450, 150, 150))) {

            assert_node(FRONT, right_front_down_cell);
            assert_node(BACK, NULL);
            assert_node(TOP, right_back_top_cell);
            assert_node(DOWN, NULL);
            assert_node(RIGHT, NULL);
            assert_node(LEFT, left_back_down_cell);
            assert_node(FRONT_TOP, right_front_top_cell);
            assert_node(FRONT_DOWN, NULL);
            assert_node(FRONT_RIGHT, NULL);
            assert_node(FRONT_LEFT, left_front_down_cell);
            assert_node(FRONT_TOP_RIGHT, NULL);
            assert_node(FRONT_TOP_LEFT, left_front_top_cell);
            assert_node(FRONT_DOWN_RIGHT, NULL);
            assert_node(FRONT_DOWN_LEFT, NULL);
            assert_node(BACK_TOP, NULL);
            assert_node(BACK_DOWN, NULL);
            assert_node(BACK_RIGHT, NULL);
            assert_node(BACK_LEFT, NULL);
            assert_node(BACK_TOP_RIGHT, NULL);
            assert_node(BACK_TOP_LEFT, NULL);
            assert_node(BACK_DOWN_RIGHT, NULL);
            assert_node(BACK_DOWN_LEFT, NULL);
            assert_node(TOP_RIGHT, NULL);
            assert_node(TOP_LEFT, left_back_top_cell);
            assert_node(DOWN_RIGHT, NULL);
            assert_node(DOWN_LEFT, NULL);
        }

        if(CENTER_EQUALS(cell->center, POINT3D(150, 450, 150))) {

            assert_node(FRONT, left_front_top_cell);
            assert_node(BACK, NULL);
            assert_node(TOP, NULL);
            assert_node(DOWN, left_back_down_cell);
            assert_node(RIGHT, right_back_top_cell);
            assert_node(LEFT, NULL);
            assert_node(FRONT_TOP, NULL);
            assert_node(FRONT_DOWN, left_front_down_cell);
            assert_node(FRONT_RIGHT, right_front_top_cell);
            assert_node(FRONT_LEFT, NULL);
            assert_node(FRONT_TOP_RIGHT, NULL);
            assert_node(FRONT_TOP_LEFT, NULL);
            assert_node(FRONT_DOWN_RIGHT, right_front_down_cell);
            assert_node(FRONT_DOWN_LEFT, NULL);
            assert_node(BACK_TOP, NULL);
            assert_node(BACK_DOWN, NULL);
            assert_node(BACK_RIGHT, NULL);
            assert_node(BACK_LEFT, NULL);
            assert_node(BACK_TOP_RIGHT, NULL);
            assert_node(BACK_TOP_LEFT, NULL);
            assert_node(BACK_DOWN_RIGHT, NULL);
            assert_node(BACK_DOWN_LEFT, NULL);
            assert_node(TOP_RIGHT, NULL);
            assert_node(TOP_LEFT, NULL);
            assert_node(DOWN_RIGHT, right_back_down_cell);
            assert_node(DOWN_LEFT, NULL);
        }

        if(CENTER_EQUALS(cell->center, POINT3D(150, 450, 450))) {

            assert_node(FRONT, NULL);
            assert_node(BACK, left_back_top_cell);
            assert_node(TOP, NULL);
            assert_node(DOWN, left_front_down_cell);
            assert_node(RIGHT, right_front_top_cell);
            assert_node(LEFT, NULL);
            assert_node(FRONT_TOP, NULL);
            assert_node(FRONT_DOWN, NULL);
            assert_node(FRONT_RIGHT, NULL);
            assert_node(FRONT_LEFT, NULL);
            assert_node(FRONT_TOP_RIGHT, NULL);
            assert_node(FRONT_TOP_LEFT, NULL);
            assert_node(FRONT_DOWN_RIGHT, NULL);
            assert_node(FRONT_DOWN_LEFT, NULL);
            assert_node(BACK_TOP, NULL);
            assert_node(BACK_DOWN, left_back_down_cell);
            assert_node(BACK_RIGHT, right_back_top_cell);
            assert_node(BACK_LEFT, NULL);
            assert_node(BACK_TOP_RIGHT, NULL);
            assert_node(BACK_TOP_LEFT, NULL);
            assert_node(BACK_DOWN_RIGHT, right_back_down_cell);
            assert_node(BACK_DOWN_LEFT, NULL);
            assert_node(TOP_RIGHT, NULL);
            assert_node(TOP_LEFT, NULL);
            assert_node(DOWN_RIGHT, right_front_down_cell);
            assert_node(DOWN_LEFT, NULL);
        }

        if(CENTER_EQUALS(cell->center, POINT3D(450, 150, 450))) {

            assert_node(FRONT, NULL);
            assert_node(BACK, right_back_down_cell);
            assert_node(TOP, right_front_top_cell);
            assert_node(DOWN, NULL);
            assert_node(RIGHT, NULL);
            assert_node(LEFT, left_front_down_cell);
            assert_node(FRONT_TOP, NULL);
            assert_node(FRONT_DOWN, NULL);
            assert_node(FRONT_RIGHT, NULL);
            assert_node(FRONT_LEFT, NULL);
            assert_node(FRONT_TOP_RIGHT, NULL);
            assert_node(FRONT_TOP_LEFT, NULL);
            assert_node(FRONT_DOWN_RIGHT, NULL);
            assert_node(FRONT_DOWN_LEFT, NULL);
            assert_node(BACK_TOP, right_back_top_cell); //
            assert_node(BACK_DOWN, NULL);
            assert_node(BACK_RIGHT, NULL);
            assert_node(BACK_LEFT, left_back_down_cell); //
            assert_node(BACK_TOP_RIGHT, NULL);
            assert_node(BACK_TOP_LEFT, left_back_top_cell); //
            assert_node(BACK_DOWN_RIGHT, NULL);
            assert_node(BACK_DOWN_LEFT, NULL);
            assert_node(TOP_RIGHT, NULL);
            assert_node(TOP_LEFT, left_front_top_cell); //
            assert_node(DOWN_RIGHT, NULL);
            assert_node(DOWN_LEFT, NULL);
        }

        if(CENTER_EQUALS(cell->center, POINT3D(450, 450, 150))) {

            assert_node(FRONT, right_front_top_cell);
            assert_node(BACK, NULL);
            assert_node(TOP, NULL);
            assert_node(DOWN, right_back_down_cell);
            assert_node(RIGHT, NULL);
            assert_node(LEFT, left_back_top_cell);
            assert_node(FRONT_TOP, NULL);
            assert_node(FRONT_DOWN, right_front_down_cell);
            assert_node(FRONT_RIGHT, NULL);
            assert_node(FRONT_LEFT, left_front_top_cell);
            assert_node(FRONT_TOP_RIGHT, NULL);
            assert_node(FRONT_TOP_LEFT, NULL);
            assert_node(FRONT_DOWN_RIGHT, NULL);
            assert_node(FRONT_DOWN_LEFT, left_front_down_cell);
            assert_node(BACK_TOP, NULL);
            assert_node(BACK_DOWN, NULL);
            assert_node(BACK_RIGHT, NULL);
            assert_node(BACK_LEFT, NULL);
            assert_node(BACK_TOP_RIGHT, NULL);
            assert_node(BACK_TOP_LEFT, NULL);
            assert_node(BACK_DOWN_RIGHT, NULL);
            assert_node(BACK_DOWN_LEFT, NULL);
            assert_node(TOP_RIGHT, NULL);
            assert_node(TOP_LEFT, NULL);
            assert_node(DOWN_RIGHT, NULL);
            assert_node(DOWN_LEFT, left_back_down_cell);
        }
    }
}

void test_custom_mesh_connectors() {

    struct grid *grid = new_grid();
    struct config *domain_config;

    domain_config = alloc_and_init_config_data();

    shput_dup_value(domain_config->config_data, "maximum_discretization", "1000.0");
    shput_dup_value(domain_config->config_data, "start_discretization", "1000.0");
    shput_dup_value(domain_config->config_data, "x_domain_limit", "128000.0");
    shput_dup_value(domain_config->config_data, "y_domain_limit", "128000.0");
    shput_dup_value(domain_config->config_data, "z_domain_limit", "128000.0");
    shput_dup_value(domain_config->config_data, "number_of_points", "202358");
    shput_dup_value(domain_config->config_data, "mesh_file", "meshes/one_ventricle.alg");

    domain_config->main_function_name = strdup("initialize_grid_with_custom_mesh");
    shput_dup_value(domain_config->config_data, "name", "Custom mesh");

    init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

    int success = ((set_spatial_domain_fn *)domain_config->main_function)(domain_config, grid);

    cr_assert(success);

    order_grid_cells(grid);

    struct config *save_mesh_config = alloc_and_init_config_data();

    save_mesh_config->main_function_name = strdup("save_as_text_or_binary");
    shput_dup_value(save_mesh_config->config_data, "output_dir", "./tests_bin");
    shput_dup_value(save_mesh_config->config_data, "print_rate", "1");
    init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_result");

    shput(save_mesh_config->config_data, "file_prefix", strdup("test_custom_mesh_connectors"));

    struct time_info ti = ZERO_TIME_INFO;

    ((save_mesh_fn *)save_mesh_config->main_function)(&ti, save_mesh_config, grid, NULL, NULL);

    free_config_data(save_mesh_config);

    struct point_3d target_cell = POINT3D(4250, 12750, 13250);
    struct point_3d target_cell2 = POINT3D(19750, 10250, 23250);
    struct point_3d target_cell3 = POINT3D(2750, 8250, 19250);

    real_cpu d = 500.0;

    FOR_EACH_CELL(grid) {

        struct cell_node *node = NULL;
        struct point_3d center = cell->center;

        if(CENTER_EQUALS(center, target_cell)) {
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, false);
            assert_node_center(RIGHT, false);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, true);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, true);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, true);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, true);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, false);
            assert_node_center(BACK_RIGHT, false);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, false);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, false);
            assert_node_center(DOWN_LEFT, true);
        }

        if(CENTER_EQUALS(center, target_cell2)) {
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, false);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, false);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, true);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, true);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, true);
            assert_node_center(BACK_TOP, false);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, false);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, false);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, false);
            assert_node_center(DOWN_LEFT, true);
        }
        if(CENTER_EQUALS(center, target_cell3)) {
            assert_node_center(FRONT, true);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, false);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, true);
            assert_node_center(FRONT_DOWN, true);
            assert_node_center(FRONT_RIGHT, true);
            assert_node_center(FRONT_LEFT, true);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, true);
            assert_node_center(FRONT_DOWN_RIGHT, true);
            assert_node_center(FRONT_DOWN_LEFT, true);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, true);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, true);
        }
    }
}

void test_cell_connector_3layer_mesh(real_cpu d) {

    struct grid *grid = new_grid();
    struct config *domain_config;

    domain_config = alloc_and_init_config_data();

    real_cpu sl = 3 * d;

    shput_dup_value(domain_config->config_data, "start_dx", sdscatprintf(sdsempty(), "%lf\n", d));
    shput_dup_value(domain_config->config_data, "start_dy", sdscatprintf(sdsempty(), "%lf\n", d));
    shput_dup_value(domain_config->config_data, "start_dz", sdscatprintf(sdsempty(), "%lf\n", d));

    domain_config->main_function_name = strdup("initialize_grid_with_cuboid_mesh");
    shput_dup_value(domain_config->config_data, "name", "Test cuboid");

    shput_dup_value(domain_config->config_data, "side_length_x", sdscatprintf(sdsempty(), "%lf\n", sl));
    shput_dup_value(domain_config->config_data, "side_length_y", sdscatprintf(sdsempty(), "%lf\n", sl));
    shput_dup_value(domain_config->config_data, "side_length_z", sdscatprintf(sdsempty(), "%lf\n", sl));

    init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

    int success = ((set_spatial_domain_fn *)domain_config->main_function)(domain_config, grid);

    cr_assert(success);

    order_grid_cells(grid);

    struct point_3d center_cell_coords = POINT3D(sl / 2.0, sl / 2.0, sl / 2.0);
    int num_cells_visited = 0;

    FOR_EACH_CELL(grid) {

        struct cell_node *node = NULL;

        struct point_3d center = cell->center;

        // Cell in the center of the mesh
        if(CENTER_EQUALS(center, center_cell_coords)) {
            num_cells_visited++;
            assert_node_center(FRONT, true);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, true);
            assert_node_center(FRONT_DOWN, true);
            assert_node_center(FRONT_RIGHT, true);
            assert_node_center(FRONT_LEFT, true);
            assert_node_center(FRONT_TOP_RIGHT, true);
            assert_node_center(FRONT_TOP_LEFT, true);
            assert_node_center(FRONT_DOWN_RIGHT, true);
            assert_node_center(FRONT_DOWN_LEFT, true);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, true);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, true);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, true);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, true);
        }

        // Going front
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, 0, 0, d))) {
            num_cells_visited++;
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, false);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, true);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, true);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, true);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, true);
        }

        // Going back
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, 0, 0, -d))) {
            num_cells_visited++;
            assert_node_center(FRONT, true);
            assert_node_center(BACK, false);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, true);
            assert_node_center(FRONT_DOWN, true);
            assert_node_center(FRONT_RIGHT, true);
            assert_node_center(FRONT_LEFT, true);
            assert_node_center(FRONT_TOP_RIGHT, true);
            assert_node_center(FRONT_TOP_LEFT, true);
            assert_node_center(FRONT_DOWN_RIGHT, true);
            assert_node_center(FRONT_DOWN_LEFT, true);
            assert_node_center(BACK_TOP, false);
            assert_node_center(BACK_DOWN, false);
            assert_node_center(BACK_RIGHT, false);
            assert_node_center(BACK_LEFT, false);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, false);
            assert_node_center(BACK_DOWN_RIGHT, false);
            assert_node_center(BACK_DOWN_LEFT, false);
            assert_node_center(TOP_RIGHT, true);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, true);
        }

        // Going top
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, 0, d, 0))) {
            num_cells_visited++;
            assert_node_center(FRONT, true);
            assert_node_center(BACK, true);
            assert_node_center(TOP, false);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, true);
            assert_node_center(FRONT_RIGHT, true);
            assert_node_center(FRONT_LEFT, true);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, true);
            assert_node_center(FRONT_DOWN_LEFT, true);
            assert_node_center(BACK_TOP, false);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, false);
            assert_node_center(BACK_DOWN_RIGHT, true);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, false);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, true);
        }

        // Going down
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, 0, -d, 0))) {
            num_cells_visited++;
            assert_node_center(FRONT, true);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, false);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, true);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, true);
            assert_node_center(FRONT_LEFT, true);
            assert_node_center(FRONT_TOP_RIGHT, true);
            assert_node_center(FRONT_TOP_LEFT, true);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, false);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, true);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, false);
            assert_node_center(BACK_DOWN_LEFT, false);
            assert_node_center(TOP_RIGHT, true);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, false);
            assert_node_center(DOWN_LEFT, false);
        }

        // Going right
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, d, 0, 0))) {
            num_cells_visited++;
            assert_node_center(FRONT, true);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, false);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, true);
            assert_node_center(FRONT_DOWN, true);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, true);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, true);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, true);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, false);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, false);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, false);
            assert_node_center(DOWN_LEFT, true);
        }

        // Going left
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, -d, 0, 0))) {
            num_cells_visited++;
            assert_node_center(FRONT, true);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, false);
            assert_node_center(FRONT_TOP, true);
            assert_node_center(FRONT_DOWN, true);
            assert_node_center(FRONT_RIGHT, true);
            assert_node_center(FRONT_LEFT, false);
            assert_node_center(FRONT_TOP_RIGHT, true);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, true);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, false);
            assert_node_center(BACK_TOP_RIGHT, true);
            assert_node_center(BACK_TOP_LEFT, false);
            assert_node_center(BACK_DOWN_RIGHT, true);
            assert_node_center(BACK_DOWN_LEFT, false);
            assert_node_center(TOP_RIGHT, true);
            assert_node_center(TOP_LEFT, false);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, false);
        }

        // Going front top
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, 0, d, d))) {
            num_cells_visited++;
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, false);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, false);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, false);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, false);
            assert_node_center(BACK_DOWN_RIGHT, true);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, false);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, true);
        }

        // Going front down
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, 0, -d, d))) {
            num_cells_visited++;
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, false);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, false);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, false);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, true);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, false);
            assert_node_center(BACK_DOWN_LEFT, false);
            assert_node_center(TOP_RIGHT, true);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, false);
            assert_node_center(DOWN_LEFT, false);
        }

        // Going front right
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, d, 0, d))) {
            num_cells_visited++;
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, false);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, false);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, false);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, true);
            assert_node_center(BACK_DOWN_RIGHT, false);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, true);
            assert_node_center(DOWN_RIGHT, false);
            assert_node_center(DOWN_LEFT, true);
        }

        // Going front left
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, -d, 0, d))) {
            num_cells_visited++;
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, true);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, false);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, false);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, true);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, false);
            assert_node_center(BACK_TOP_RIGHT, true);
            assert_node_center(BACK_TOP_LEFT, false);
            assert_node_center(BACK_DOWN_RIGHT, true);
            assert_node_center(BACK_DOWN_LEFT, false);
            assert_node_center(TOP_RIGHT, true);
            assert_node_center(TOP_LEFT, false);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, false);
        }

        // Going front top right
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, d, d, d))) {
            num_cells_visited++;
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, false);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, false);
            assert_node_center(LEFT, true);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, false);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, false);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, false);
            assert_node_center(BACK_LEFT, true);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, false);
            assert_node_center(BACK_DOWN_RIGHT, false);
            assert_node_center(BACK_DOWN_LEFT, true);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, false);
            assert_node_center(DOWN_RIGHT, false);
            assert_node_center(DOWN_LEFT, true);
        }

        // Going front top left
        if(CENTER_EQUALS(center, TRANSLATE(center_cell_coords, -d, d, d))) {
            num_cells_visited++;
            assert_node_center(FRONT, false);
            assert_node_center(BACK, true);
            assert_node_center(TOP, false);
            assert_node_center(DOWN, true);
            assert_node_center(RIGHT, true);
            assert_node_center(LEFT, false);
            assert_node_center(FRONT_TOP, false);
            assert_node_center(FRONT_DOWN, false);
            assert_node_center(FRONT_RIGHT, false);
            assert_node_center(FRONT_LEFT, false);
            assert_node_center(FRONT_TOP_RIGHT, false);
            assert_node_center(FRONT_TOP_LEFT, false);
            assert_node_center(FRONT_DOWN_RIGHT, false);
            assert_node_center(FRONT_DOWN_LEFT, false);
            assert_node_center(BACK_TOP, false);
            assert_node_center(BACK_DOWN, true);
            assert_node_center(BACK_RIGHT, true);
            assert_node_center(BACK_LEFT, false);
            assert_node_center(BACK_TOP_RIGHT, false);
            assert_node_center(BACK_TOP_LEFT, false);
            assert_node_center(BACK_DOWN_RIGHT, true);
            assert_node_center(BACK_DOWN_LEFT, false);
            assert_node_center(TOP_RIGHT, false);
            assert_node_center(TOP_LEFT, false);
            assert_node_center(DOWN_RIGHT, true);
            assert_node_center(DOWN_LEFT, false);
        }
    }

    cr_assert_eq(13, num_cells_visited);

    struct config *save_mesh_config = alloc_and_init_config_data();

    save_mesh_config->main_function_name = strdup("save_as_text_or_binary");
    shput_dup_value(save_mesh_config->config_data, "output_dir", "./tests_bin");
    shput_dup_value(save_mesh_config->config_data, "print_rate", "1");
    init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_result");

    sds file_prefix = sdscatprintf(sdsempty(), "test_connectors_%lf_%lf", d, sl);

    shput(save_mesh_config->config_data, "file_prefix", strdup(file_prefix));

    struct time_info ti = ZERO_TIME_INFO;

    ((save_mesh_fn *)save_mesh_config->main_function)(&ti, save_mesh_config, grid, NULL, NULL);

    free_config_data(save_mesh_config);
}

Test(cell_conectors, cuboid_26_neighbours_50um) {
    test_cell_connector_3layer_mesh(50.0);
}

Test(cell_conectors, cuboid_26_neighbours_100um) {
    test_cell_connector_3layer_mesh(100.0);
}

Test(cell_conectors, cuboid_26_neighbours_200um) {
    test_cell_connector_3layer_mesh(200.0);
}

Test(cell_conectors, cuboid_26_neighbours_300um) {
    test_cell_connector_3layer_mesh(300.0);
}

Test(cell_conectors, custom_mesh) {
    test_custom_mesh_connectors();
}
