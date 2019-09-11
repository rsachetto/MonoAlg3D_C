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

    // TODO: Consider also the diameter of the Purkinje cell ...
    real_cpu side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu,side_length, config->config_data, "start_discretization");
    
    char *network_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(network_file,config->config_data,"network_file");

    char *celular_model_name = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(celular_model_name,config->config_data,"celular_model");

    init_purkinje_ode_solver_with_cell_model(the_ode_solver,celular_model_name);

    // TODO: Implement this test
    int success = check_purkinje_input();

    if (!success)
    {
        return 0;
    }

    print_to_stdout_and_file("Loading a custom Purkinje Network:> %s\n", name);
    print_to_stdout_and_file("Using the Purkinje library function:> \"initialize_purkinje_with_custom_mesh\"\n");
    print_to_stdout_and_file("Discretization for the Purkinje Network Mesh:> %g um\n",side_length);
    print_to_stdout_and_file("Celular model for the Purkinje :> %s\n",the_ode_solver->model_data.model_library_path);
    set_custom_purkinje_network(the_grid->the_purkinje, network_file, side_length);

    // Populate the 'purkinje_cells' linked-list with the nodes from the graph
    //        Some parameters from the 'cell_node' structure will not be used
    initialize_and_construct_grid_purkinje(the_grid);

    free (network_file);

    return 1;
}

// TO DO: Build some benchmark Purkinje network for quick tests ...