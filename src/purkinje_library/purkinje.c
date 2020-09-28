//
// Created by bergolho on 19/07/18.
//

#include "purkinje_helpers.h"

#include "../config/purkinje_config.h"
#include "../monodomain/constants.h"
#include "../libraries_common/common_data_structures.h"
#include "../config_helpers/config_helpers.h"
#include "../3dparty/sds/sds.h"
#include "../utils/file_utils.h"
#include "../3dparty/stb_ds.h"
#include "../utils/utils.h"

#include <assert.h>
#include <time.h>
#include <unistd.h>

SET_SPATIAL_PURKINJE (initialize_purkinje_with_custom_mesh)  {

    printf("On 'initialize_purkinje_with_custom_mesh'\n");

    char *name = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(name, config->config_data, "name");

    real_cpu dx = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu,dx, config->config_data, "dx");
    
    char *network_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(network_file,config->config_data,"network_file");

    log_to_stdout_and_file("Loading a custom Purkinje Network:> %s\n", name);
    log_to_stdout_and_file("Using the Purkinje library function:> \"initialize_purkinje_with_custom_mesh\"\n");
    log_to_stdout_and_file("Discretization for the Purkinje Network Mesh:> %g um\n",dx);
    log_to_stdout_and_file("Celular model for the Purkinje :> %s\n",the_ode_solver->model_data.model_library_path);
    set_custom_purkinje_network(the_grid->purkinje, network_file, dx);

    // Populate the 'purkinje_cells' linked-list with the nodes from the graph
    //        Some parameters from the 'cell_node' structure will not be used
    initialize_and_construct_grid_purkinje(the_grid);

    free (network_file);

    // Before returning test if there is an error in Purkinje mesh
    int success = check_purkinje_mesh_for_errors(the_grid->purkinje->network);
    if (!success)
    {
        return 0;
    }

    return 1;
}
