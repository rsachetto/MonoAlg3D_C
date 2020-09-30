
// TODO: Include all the Purkinje related functions to save meshes here ...
//       New library --> 'libdefault_save_mesh_purkinje.so'

/*
SAVE_MESH(save_as_vtp_purkinje) {

    char *output_dir;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");

    if(!initialized) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config->config_data, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config->config_data, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config->config_data, "compression_level");

        if(compress) binary = true;

        initialized = true;
    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[2], config->config_data, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[4], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[5], config->config_data, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[0], config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[1], config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[2], config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[3], config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[4], config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[5], config->config_data, "max_z");
    }

    write_transmembrane_potential_vtp(&vtk_polydata, the_grid, output_dir, file_prefix, time_info->iteration, time_info->current_t,\
                                             save_pvd, compress, compression_level, binary, clip_with_plain, plain_coords, clip_with_bounds, bounds);
}
*/