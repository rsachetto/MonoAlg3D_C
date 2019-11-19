//
// Created by bergolho on 19/07/18.
//

#include "purkinje_helpers.h"

#include "../config/purkinje_config.h"
#include "../monodomain/constants.h"
#include "../libraries_common/common_data_structures.h"
#include "../config_helpers/config_helpers.h"
#include "../string/sds.h"
#include "../utils/file_utils.h"
#include "../single_file_libraries/stb_ds.h"
#include "../utils/utils.h"

#include <assert.h>
#include <time.h>
#include <unistd.h>

SET_SPATIAL_PURKINJE (initialize_purkinje_with_custom_mesh)
{
    printf("On 'initialize_purkinje_with_custom_mesh'\n");

    char *name = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(name, config->config_data, "name");

    real_cpu side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu,side_length, config->config_data, "start_discretization");

    char *network_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(network_file,config->config_data,"network_file");

    real_cpu rpmj = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu,rpmj, config->config_data, "rpmj");

    real_cpu pmj_scale = 0.1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu,pmj_scale, config->config_data, "pmj_scale");

    bool calc_retro_propagation = true;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(bool,calc_retro_propagation, config->config_data, "retro_propagation");

    init_ode_solver_with_cell_model(the_ode_solver);

    // TODO: Implement this test
    int success = check_purkinje_input();

    if (!success)
    {
        return 0;
    }

    print_to_stdout_and_file("Loading a custom Purkinje Network:> %s\n", name);
    print_to_stdout_and_file("Using the Purkinje library function:> \"initialize_purkinje_with_custom_mesh\"\n");
    print_to_stdout_and_file("Discretization for the Purkinje Network Mesh:> %g um\n",side_length);
    print_to_stdout_and_file("Purkinje-Muscle-Junction resistance:> %g um\n",rpmj);
    print_to_stdout_and_file("Purkinje-Muscle-Junction scale:> %g\n",pmj_scale);
    print_to_stdout_and_file("Celular model for the Purkinje :> %s\n",the_ode_solver->model_data.model_library_path);

    print_to_stdout_and_file("Loading Purkinje mesh ...\n");
    set_custom_purkinje_network(the_grid->the_purkinje, network_file, side_length, rpmj, pmj_scale, calc_retro_propagation);

    // Populate the 'purkinje_cells' linked-list with the nodes from the graph
    //        Some parameters from the 'cell_node' structure will not be used
    initialize_and_construct_grid_purkinje(the_grid);

    free (network_file);

    return 1;
}

// TO DO: Build some benchmark Purkinje network for quick tests ...
