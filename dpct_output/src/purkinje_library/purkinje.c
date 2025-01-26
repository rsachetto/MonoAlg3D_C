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
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(name, config, "name");

    real_cpu dx = 100.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu,dx, config, "dx");

    char *network_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(network_file,config,"network_file");

    log_info("Loading a custom Purkinje Network:> %s\n", name);
    log_info("Using the Purkinje library function:> \"initialize_purkinje_with_custom_mesh\"\n");
    log_info("Discretization for the Purkinje Network Mesh:> %g um\n",dx);
    log_info("Celular model for the Purkinje :> %s\n",the_ode_solver->model_data.model_library_path);
    set_custom_purkinje_network(the_grid->purkinje, network_file, dx);

    // Populate the 'purkinje_cells' linked-list with the nodes from the graph
    //        Some parameters from the 'cell_node' structure will not be used
    initialize_and_construct_grid_purkinje(the_grid);

    free (network_file);

    // Before returning test if there is an error in Purkinje mesh
    int success = check_purkinje_mesh_for_errors(the_grid->purkinje->network);

    if (!success) {
        return 0;
    }

    return 1;
}

SET_SPATIAL_PURKINJE (initialize_purkinje_coupling_with_custom_mesh) {

    char *name = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(name, config, "name");

    real_cpu dx = 100.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, dx, config, "dx");

    real_cpu rpmj = 1000.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, rpmj, config, "rpmj");

    real_cpu pmj_scale = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, pmj_scale, config, "pmj_scale");

    real_cpu asymm_ratio = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, asymm_ratio, config, "asymm_ratio");

    uint32_t nmin_pmj = 10;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, nmin_pmj, config, "nmin_pmj");

    uint32_t nmax_pmj = 30;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, nmax_pmj, config, "nmax_pmj");

    char *network_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(network_file,config,"network_file");

    bool retro_propagation = true;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(retro_propagation, config, "retro_propagation");

    char *pmj_location_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(pmj_location_file, config, "pmj_location_file");

    log_info("Loading a custom Purkinje Network:> %s\n", name);
    log_info("Using the Purkinje library function:> \"initialize_purkinje_with_custom_mesh\"\n");
    log_info("Discretization for the Purkinje Network Mesh:> %g um\n",dx);
    log_info("Purkinje-Muscle-Junction resistance:> %g kohm\n",rpmj);
    log_info("Purkinje-Muscle-Junction scale:> %g\n",pmj_scale);
    log_info("Purkinje-Muscle-Junction asymmetry conduction ratio:> %g\n",asymm_ratio);
    log_info("Minimum tissue nodes inside Purkinje-Muscle-Junction:> %u\n",nmin_pmj);
    log_info("Maximum tissue nodes inside Purkinje-Muscle-Junction:> %u\n",nmax_pmj);
    log_info("Cellular model for the Purkinje :> %s\n",the_ode_solver->model_data.model_library_path);
    set_custom_purkinje_network(the_grid->purkinje, network_file, dx);

    set_purkinje_coupling_parameters(the_grid->purkinje->network,rpmj,pmj_scale,asymm_ratio,nmin_pmj,nmax_pmj,retro_propagation,pmj_location_file);

    // Populate the 'purkinje_cells' linked-list with the nodes from the graph
    //        Some parameters from the 'cell_node' structure will not be used
    initialize_and_construct_grid_purkinje(the_grid);

    free (network_file);

    // Before returning test if there is an error in Purkinje mesh
    int success = check_purkinje_mesh_for_errors(the_grid->purkinje->network);

    if (!success) {
        return 0;
    }

    return 1;
}
