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

SET_SPATIAL_PURKINJE (initialize_purkinje_with_custom_mesh) {

    char *name = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(name, config->config_data, "name");

    real_cpu side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length, config->config_data, "start_discretization");
    
    char *network_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(network_file,config->config_data,"network_file");

    log_to_stdout_and_file("Loading a custom Purkinje Network:> %s\n", name);
    log_to_stdout_and_file("Using the Purkinje library function:> \"initialize_purkinje_with_custom_mesh\"\n");
    log_to_stdout_and_file("Discretization for the Purkinje Network Mesh:> %g um\n",side_length);
    log_to_stdout_and_file("Cellular model for the Purkinje :> %s\n",the_ode_solver->model_data.model_library_path);
    set_custom_purkinje_network(the_grid->purkinje, network_file, side_length);

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

SET_SPATIAL_PURKINJE (initialize_purkinje_coupling_with_custom_mesh) {

    char *name = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(name, config->config_data, "name");

    real_cpu side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length, config->config_data, "start_discretization");
    
    real_cpu rpmj = 1000.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, rpmj, config->config_data, "rpmj");

    real_cpu pmj_scale = 1000.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, pmj_scale, config->config_data, "pmj_scale");

    real_cpu asymm_ratio = 0.01;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, asymm_ratio, config->config_data, "asymm_ratio");

    uint32_t nmin_pmj = 10;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, nmin_pmj, config->config_data, "nmin_pmj");

    uint32_t nmax_pmj = 30;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, nmax_pmj, config->config_data, "nmax_pmj");

    char *network_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(network_file,config->config_data,"network_file");

    bool retro_propagation = true;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(retro_propagation, config->config_data, "retro_propagation");

    char *pmj_location_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(pmj_location_file, config->config_data, "pmj_location_file");

    log_to_stdout_and_file("Loading a custom Purkinje Network:> %s\n", name);
    log_to_stdout_and_file("Using the Purkinje library function:> \"initialize_purkinje_with_custom_mesh\"\n");
    log_to_stdout_and_file("Discretization for the Purkinje Network Mesh:> %g um\n",side_length);
    log_to_stdout_and_file("Purkinje-Muscle-Junction resistance:> %g kohm\n",rpmj);
    log_to_stdout_and_file("Purkinje-Muscle-Junction scale:> %g\n",pmj_scale);
    log_to_stdout_and_file("Purkinje-Muscle-Junction asymmetry conduction ratio:> %g\n",asymm_ratio);
    log_to_stdout_and_file("Minimum tissue nodes inside Purkinje-Muscle-Junction:> %g\n",nmin_pmj);
    log_to_stdout_and_file("Maximum tissue nodes inside Purkinje-Muscle-Junction:> %g\n",nmax_pmj);
    log_to_stdout_and_file("Cellular model for the Purkinje :> %s\n",the_ode_solver->model_data.model_library_path);
    set_custom_purkinje_network(the_grid->purkinje, network_file, side_length);

    set_purkinje_coupling_parameters(the_grid->purkinje->network,rpmj,pmj_scale,asymm_ratio,nmin_pmj,nmax_pmj,retro_propagation,pmj_location_file);

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