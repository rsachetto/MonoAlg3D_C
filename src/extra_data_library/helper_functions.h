//
// Created by sachetto on 19/03/2021.
//

#ifndef MONOALG3D_C_EXTRA_DATA_HELPER_FUNCTIONS_H
#define MONOALG3D_C_EXTRA_DATA_HELPER_FUNCTIONS_H

#include <unistd.h>
#include "../common_types/common_types.h"
#include "../config/config_common.h"

#define SET_EXTRA_DATA_SIZE(value) *extra_data_size = (value)

struct extra_data_for_fibrosis {
    real atpi;
    real Ko;
    real Ki;
    real GNa_multiplicator;
    real GCaL_multiplicator;
    real INaCa_multiplicator;
    real Vm_modifier;
    real *fibrosis;
};

struct extra_data_for_torord {
    real INa_Multiplier; 
    real ICaL_Multiplier;
    real Ito_Multiplier;
    real INaL_Multiplier;
    real IKr_Multiplier; 
    real IKs_Multiplier; 
    real IK1_Multiplier; 
    real IKb_Multiplier; 
    real INaCa_Multiplier;
    real INaK_Multiplier;  
    real INab_Multiplier;  
    real ICab_Multiplier;  
    real IpCa_Multiplier;  
    real ICaCl_Multiplier;
    real IClb_Multiplier; 
    real Jrel_Multiplier; 
    real Jup_Multiplier;
    real *initial_ss_endo;
    real *initial_ss_epi;
    real *initial_ss_mid;
    real *transmurality;
};

struct extra_data_for_trovato {
    real GNa_Multiplier;
    real GNaL_Multiplier;
    real GCaT_Multiplier;
    real Gto_Multiplier;
    real Gsus_Multiplier;
    real Gkr_Multiplier;
    real Gks_Multiplier;
    real GfNa_Multiplier;
    real GfK_Multiplier;
    real GK1_Multiplier;
    real GNCX_Multiplier;
    real GNaK_Multiplier;
    real INa_Multiplier; 
    real ICaL_Multiplier;
    real ICaNa_Multiplier;
    real ICaK_Multiplier;
    real Ito_Multiplier;
    real INaL_Multiplier;
    real IKr_Multiplier; 
    real IKs_Multiplier; 
    real IK1_Multiplier; 
    real INaCa_Multiplier;
    real INaK_Multiplier;  
    real INab_Multiplier;  
    real ICab_Multiplier;
    real ICaT_Multiplier;
    real Isus_Multiplier;
    real If_Multiplier;
    real IpCa_Multiplier;
    //real *purkinje_tags;
};

struct extra_data_for_fibrosis * set_common_schemia_data(struct config *config, uint32_t num_cells);
struct extra_data_for_torord * set_common_torord_data (struct config *config, uint32_t num_cells);
struct extra_data_for_torord * set_common_torord_dyncl_data (struct config *config, uint32_t num_cells);
struct extra_data_for_trovato * set_common_trovato_data (struct config *config, uint32_t num_cells);

#endif // MONOALG3D_C_EXTRA_DATA_HELPER_FUNCTIONS_H
