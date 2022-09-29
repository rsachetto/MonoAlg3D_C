#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "3dparty/sds/sds.h"
#include "3dparty/stb_ds.h"
#include "3dparty/ini_parser/ini.h"
#include "config/config_parser.h"
#include "utils/file_utils.h"
#include "vtk_utils/vtk_unstructured_grid.h"

static void convert_file(const char *input, const char *output, const char *file_name, struct string_voidp_hash_entry *extra_data_config) {

    sds full_input_path;

    static struct vtk_unstructured_grid *vtk_grid = NULL;
    sds full_output_path;

    full_input_path = sdsnew(input);

    struct path_information input_file_info;

    if(file_name != NULL) {
        if(!ENDS_WITH_SLASH(full_input_path)) {
            full_input_path = sdscat(full_input_path, "/");
        }

        full_input_path = sdscat(full_input_path, file_name);
    }

    get_path_information(full_input_path, &input_file_info);

    static int count = 0;

    if(FILE_HAS_EXTENSION(input_file_info, "txt") || FILE_HAS_EXTENSION(input_file_info, "alg") ||
       FILE_HAS_EXTENSION(input_file_info, "geo") || FILE_HAS_EXTENSION_PREFIX(input_file_info, "Esca")) {

        if(vtk_grid == NULL) {
            vtk_grid = new_vtk_unstructured_grid_from_file(full_input_path, false);
            if(FILE_HAS_EXTENSION(input_file_info, "geo")) {
                return;
            }
        }

        if(!vtk_grid) {
            fprintf(stderr, "%s is not a valid simulation file. Skipping!!\n", full_input_path);
        } else {

            full_output_path = sdsnew(output);

            if(!ENDS_WITH_SLASH(full_output_path)) {
                full_output_path = sdscat(full_output_path, "/");
            }

            if(FILE_HAS_EXTENSION_PREFIX(input_file_info, "Esca")) {
                set_vtk_grid_values_from_ensight_file(vtk_grid, full_input_path);
            }

            char ext[4] = ".vtu";
            bool alg = false;

            if(FILE_HAS_EXTENSION(input_file_info, "alg")) {
                ext[1] = 'v';
                ext[2] = 't';
                ext[3] = 'k';
                alg = true;
            }

            if(FILE_HAS_EXTENSION_PREFIX(input_file_info, "Esca")) {
                full_output_path = sdscatfmt(full_output_path, "V_it_%i%s", count, ext);
                count++;
            } else {
                full_output_path = sdscat(full_output_path, input_file_info.filename_without_extension);
                full_output_path = sdscat(full_output_path, ext);
            }

            printf("Converting %s to %s\n", full_input_path, full_output_path);

            if(alg){
                save_vtk_unstructured_grid_as_legacy_vtk(vtk_grid, full_output_path, false, false, extra_data_config);
            } else {
                save_vtk_unstructured_grid_as_vtu_compressed(vtk_grid, full_output_path, 6);
            }

            if(FILE_HAS_EXTENSION(input_file_info, "txt") || FILE_HAS_EXTENSION(input_file_info, "alg")) {
                free_vtk_unstructured_grid(vtk_grid);
                vtk_grid = NULL;
            }

            sdsfree(full_output_path);

            free_path_information(&input_file_info);
            sdsfree(full_input_path);
        }

    } else if(FILE_HAS_EXTENSION(input_file_info, "vtu")) {

        vtk_grid = new_vtk_unstructured_grid_from_file(full_input_path, false);

        if(!vtk_grid) {
            fprintf(stderr, "%s is not a valid simulation file. Skipping!!\n", full_input_path);
        } else {
            full_output_path = sdsnew(output);

            if(!ENDS_WITH_SLASH(full_output_path)) {
                full_output_path = sdscat(full_output_path, "/");
            }

            full_output_path = sdscat(full_output_path, input_file_info.filename_without_extension);
            full_output_path = sdscat(full_output_path, ".txt");

            printf("Converting %s to %s\n", full_input_path, full_output_path);

            save_vtk_unstructured_grid_as_alg_file(vtk_grid, full_output_path, false);
            free_vtk_unstructured_grid(vtk_grid);
            sdsfree(full_output_path);

            free_path_information(&input_file_info);
            sdsfree(full_input_path);
        }
    }
}

#define SET_OUT_DIR(input)                                                                                                                                     \
    if(ENDS_WITH_SLASH(input)) {                                                                                                                               \
        output = sdscatfmt(output, "%sconverted_files", input);                                                                                                \
    } else {                                                                                                                                                   \
        output = sdscatfmt(output, "%s/converted_files", input);                                                                                               \
    }                                                                                                                                                          \
    create_dir(output);

int main(int argc, char **argv) {

    struct conversion_options *options = new_conversion_options();

    parse_conversion_options(argc, argv, options);

    char *input = options->input;
    char *output = options->output;

    struct path_information input_info;

    get_path_information(input, &input_info);

    if(options->conversion_config_file) {
        if (ini_parse(options->conversion_config_file, parse_converter_config_file, options) < 0) {
            fprintf(stderr, "Error parsing config file %s\n", options->conversion_config_file);
        }
    }

    if(!input_info.exists) {
        fprintf(stderr,
                "Invalid directory or file (%s)! The input parameter should be and valid simulation file, alg file or a directory containing simulations!\n",
                input);
        return EXIT_FAILURE;
    }

    if(!output) {

        output = sdsempty();

        if(input_info.is_dir) {

            SET_OUT_DIR(input);

            string_array geo_file = list_files_from_dir(input, NULL, "geo", NULL, true);
            string_array files_list = NULL;
            if(arrlen(geo_file) > 0) {

                convert_file(input, output, geo_file[0], options->extra_data_config);
                files_list = list_files_from_dir(input, "Vm.", NULL, NULL, true);
                int num_files = arrlen(files_list);

                if(num_files == 0) {
                    fprintf(stderr, "Directory %s is empty\n", input);
                    exit(EXIT_FAILURE);
                } else {
                    for(int i = 0; i < num_files; i++) {
                        convert_file(input, output, files_list[i], options->extra_data_config);
                    }
                }

            } else {

                files_list = list_files_from_dir(input, NULL, NULL, NULL, true);
                int num_files = arrlen(files_list);

                if(num_files == 0) {
                    fprintf(stderr, "Directory %s is empty\n", input);
                    exit(EXIT_FAILURE);
                } else {
                    for(int i = 0; i < num_files; i++) {
                        convert_file(input, output, files_list[i], options->extra_data_config);
                    }
                }
            }
        } else {

            SET_OUT_DIR(input_info.dir_name);

            create_dir(output);
            convert_file(input, output, NULL, options->extra_data_config);
        }
    }

    return EXIT_SUCCESS;
}
