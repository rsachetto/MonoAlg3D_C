#include "utils/file_utils.h"

#include "config/config_parser.h"
#include "vtk_utils/vtk_unstructured_grid.h"

#include "3dparty/sds/sds.h"
#include "3dparty/stb_ds.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/stat.h>

int main(int argc, char **argv) {

    struct conversion_options *options = new_conversion_options();

    parse_conversion_options(argc, argv, options);

    char *input = options->input;
    char *output = options->output;

    struct path_information input_info;

    get_path_information(input, &input_info);

    if(!input_info.exists) {
        fprintf(stderr, "Invalid directory or txt file (%s)! The input parameter should be and valid txt simulation file or a directory containing simulations!\n", input);
    }
    else if(input_info.is_dir) {

        if(!output) {

            output = sdsempty();

            if(ENDS_WITH_SLASH(input)) {
                output = sdscatfmt(output, "%sconverted_files", input);
            } else {
                output = sdscatfmt(output, "%s/converted_files", input);
            }
        }

        create_dir(output);

        string_array files_list = list_files_from_dir_sorted(input, NULL);

        int num_files = arrlen(files_list);

        if(num_files == 0) {
            fprintf(stderr, "Directory %s is empty\n", input);
            exit(EXIT_FAILURE);
        }
        else {

            struct path_information file_info;
            sds full_input_path;
            sds full_output_path;
            struct vtk_unstructured_grid *vtk_grid = NULL;

            for(int i = 0; i < num_files; i++) {

                full_input_path = sdsnew(input);

                if(!ENDS_WITH_SLASH(full_input_path)) {
                    full_input_path = sdscat(full_input_path, "/");
                }

                full_input_path = sdscat(full_input_path, files_list[i]);

                get_path_information(full_input_path, &file_info);

                if(FILE_HAS_EXTENSION(file_info.file_extension, "txt")) {

                    //TODO: maybe this should be a function
                    vtk_grid = new_vtk_unstructured_grid_from_file(full_input_path);

                    if(!vtk_grid) {
                        fprintf(stderr, "%s is not a valid simulation file. Skipping!!\n", full_input_path);
                    }
                    else {

                        full_output_path = sdsnew(output);

                        if(!ENDS_WITH_SLASH(full_output_path)) {
                            full_output_path = sdscat(full_output_path, "/");
                        }

                        full_output_path = sdscat(full_output_path, file_info.filename_without_extension);
                        full_output_path = sdscat(full_output_path, ".vtu");

                        printf("Converting %s to %s\n", full_input_path, full_output_path);

                        save_vtk_unstructured_grid_as_vtu_compressed(vtk_grid, full_output_path, 6);
                        free_vtk_unstructured_grid(vtk_grid);
                        sdsfree(full_output_path);
                    }

                }
                else if(FILE_HAS_EXTENSION(file_info.file_extension, "vtu")) {

                    vtk_grid = new_vtk_unstructured_grid_from_file(full_input_path);

                    if(!vtk_grid) {
                        fprintf(stderr, "%s is not a valid simulation file. Skipping!!\n", full_input_path);
                    }
                    else {
                        full_output_path = sdsnew(output);

                        if(!ENDS_WITH_SLASH(full_output_path)) {
                            full_output_path = sdscat(full_output_path, "/");
                        }

                        full_output_path = sdscat(full_output_path, file_info.filename_without_extension);
                        full_output_path = sdscat(full_output_path, ".txt");

                        printf("Converting %s to %s\n", full_input_path, full_output_path);

                        save_vtk_unstructured_grid_as_alg_file(vtk_grid, full_output_path, false);
                        free_vtk_unstructured_grid(vtk_grid);
                        sdsfree(full_output_path);
                    }
                }


                free_path_information(&file_info);
                sdsfree(full_input_path);
            }

        }
    }


    return EXIT_SUCCESS;
}
