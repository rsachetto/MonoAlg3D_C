////
//// Created by sachetto on 06/10/17.
////
#include <criterion/criterion.h>
#include <signal.h>

#include "../alg/grid/grid.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../config/config_parser.h"
#include "../utils/file_utils.h"
#include "../3dparty/ini_parser/ini.h"
#include "../3dparty/sds/sds.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

int test_cuboid_mesh(char *start_dx, char* start_dy, char* start_dz, char* side_length_x, char* side_length_y, char* side_length_z, bool save, bool compress,  bool binary, int id) {

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

    int success = ((set_spatial_domain_fn*)domain_config->main_function)(domain_config, grid);

    if(!success ) {
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
            if (cell->center.x > max_x) {
                max_x = cell->center.x;
            }

            if (cell->center.y > max_y) {
                max_y = cell->center.y;
            }

            if (cell->center.z > max_z) {
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
        }
        else if(binary) {
            shput(save_mesh_config->config_data, "binary", strdup("yes"));
        }

        shput(save_mesh_config->config_data, "save_pvd", strdup("no"));

        struct time_info ti = ZERO_TIME_INFO;

        ((init_save_mesh_fn *)save_mesh_config->init_function)(save_mesh_config);
        ((save_mesh_fn *)save_mesh_config->main_function)(&ti, save_mesh_config, grid, NULL);
        ((end_save_mesh_fn *)save_mesh_config->end_function)(save_mesh_config);

        free_config_data(save_mesh_config);
        sdsfree(file_prefix);

    }

    cr_assert_float_eq(max_x+(start_dx_value/2.0), atof(side_length_x), 1e-16);
    cr_assert_float_eq(max_y+(start_dy_value/2.0), atof(side_length_y), 1e-16);
    cr_assert_float_eq(max_z+(start_dz_value/2.0), atof(side_length_z), 1e-16);
    cr_assert_eq(nx*ny*nz, grid->num_active_cells);

    clean_and_free_grid(grid);
    free_config_data(domain_config);

    return 1;

}

int compare_two_binary_files(FILE *fp1, FILE *fp2)
{
    char ch1, ch2;
    int flag = 0;

    while (((ch1 = fgetc(fp1)) != EOF) &&((ch2 = fgetc(fp2)) != EOF))
    {
        /*
          * character by character comparision
          * if equal then continue by comparing till the end of files
          */
        if (ch1 == ch2)
        {
            flag = 1;
            continue;
        }
            /*
              * If not equal then returns the byte position
              */
        else
        {
            fseek(fp1, -1, SEEK_CUR);
            flag = 0;
            break;
        }
    }

    if (flag == 0)
    {
        return ftell(fp1)+1;
    }
    else
    {
        return -1;
    }
}

Test (mesh_load, cuboid_mesh_100_100_100_1000_1000_1000) {
    int success  = test_cuboid_mesh("100", "100", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_200_100_100_1000_1000_1000) {
    int success = test_cuboid_mesh("200", "100", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_100_200_100_1000_1000_1000) {

    int success = test_cuboid_mesh("100", "200", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_100_100_200_1000_1000_1000) {
    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);

}

Test (mesh_load, cuboid_mesh_100_100_100_1000_1000_2000) {
    int success = test_cuboid_mesh("100", "100", "100", "1000", "1000", "2000", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_150_150_150_1500_1500_1500) {
    int success = test_cuboid_mesh("150", "150", "150", "1500", "1500", "1500", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_150_150_150_1500_1500_3000) {
    int success  = test_cuboid_mesh("150", "150", "150", "1500", "1500", "3000", false, false, false, 0);
    cr_assert(success);
}
Test (mesh_load, cuboid_mesh_300_150_150_1500_1500_3000) {
    int success = test_cuboid_mesh("300", "150", "150", "1500", "1500", "3000", false, false, false, 0);
    cr_assert(!success);
}

Test (mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_binary) {

    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", true, false, true, 1);
    cr_assert(success);

    FILE *f1 = fopen("tests_bin/test_100_100_200_1000_1000_1000_1_it_0.vtu", "r");
    FILE *f2 = fopen("tests_bin/gold_vtu_mesh_binary.vtu", "r");

    success = compare_two_binary_files(f1, f2);

    fclose(f1);
    fclose(f2);

    cr_assert(success == -1);
}

Test (mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_compressed) {

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

Test (mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_plain) {

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