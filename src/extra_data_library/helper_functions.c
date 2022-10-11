//
// Created by sachetto on 19/03/2021.
//

#include "helper_functions.h"
#include <stdlib.h>
#include "../config_helpers/config_helpers.h"


struct extra_data_for_fibrosis* set_common_schemia_data(struct config *config, uint32_t num_cells) {

    struct extra_data_for_fibrosis *extra_data = MALLOC_ONE_TYPE(struct extra_data_for_fibrosis);

    real atpi = 6.8;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, atpi, config, "atpi");

    real Ko = 5.4;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ko, config, "Ko");

    real Ki = 138.3;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ki, config, "Ki");

    real GNa_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GNa_multiplicator, config, "GNa_multiplicator");

    real GCaL_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GCaL_multiplicator, config, "GCaL_multiplicator");

    real INaCa_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaCa_multiplicator, config, "INaCa_multiplicator");

    real Vm_modifier = 0.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Vm_modifier, config, "Vm_modifier");

    extra_data->atpi = atpi;
    extra_data->Ko = Ko;
    extra_data->Ki = Ki;
    extra_data->GNa_multiplicator = GNa_multiplicator;
    extra_data->GCaL_multiplicator = GCaL_multiplicator;
    extra_data->INaCa_multiplicator = INaCa_multiplicator;
    extra_data->Vm_modifier = Vm_modifier;
    extra_data->fibrosis = MALLOC_ARRAY_OF_TYPE(real, num_cells);

    return extra_data;
}

struct extra_data_for_torord * set_common_torord_data (struct config *config, uint32_t num_cells) {
    struct extra_data_for_torord *extra_data = MALLOC_ONE_TYPE(struct extra_data_for_torord);

    real INa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INa_Multiplier, config, "INa_Multiplier");

    real ICaL_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaL_Multiplier, config, "ICaL_Multiplier");

    real Ito_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ito_Multiplier, config, "Ito_Multiplier");

    real INaL_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaL_Multiplier, config, "INaL_Multiplier");

    real IKr_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKr_Multiplier, config, "IKr_Multiplier");

    real IKs_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKs_Multiplier, config, "IKs_Multiplier");

    real IK1_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IK1_Multiplier, config, "IK1_Multiplier");

    real IKb_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKb_Multiplier, config, "IKb_Multiplier");

    real INaCa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaCa_Multiplier, config, "INaCa_Multiplier");

    real INaK_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaK_Multiplier, config, "INaK_Multiplier");

    real INab_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INab_Multiplier, config, "INab_Multiplier");

    real ICab_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICab_Multiplier, config, "ICab_Multiplier");

    real IpCa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IpCa_Multiplier, config, "IpCa_Multiplier");

    real ICaCl_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaCl_Multiplier, config, "ICaCl_Multiplier");

    real IClb_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IClb_Multiplier, config, "IClb_Multiplier");

    real Jrel_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Jrel_Multiplier, config, "Jrel_Multiplier");

    real Jup_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Jup_Multiplier, config, "Jup_Multiplier");

    extra_data->INa_Multiplier = INa_Multiplier;
    extra_data->ICaL_Multiplier = ICaL_Multiplier;
    extra_data->Ito_Multiplier = Ito_Multiplier;
    extra_data->INaL_Multiplier = INaL_Multiplier;
    extra_data->IKr_Multiplier = IKr_Multiplier;
    extra_data->IKs_Multiplier = IKs_Multiplier;
    extra_data->IK1_Multiplier = IK1_Multiplier;
    extra_data->IKb_Multiplier = IKb_Multiplier;
    extra_data->INaCa_Multiplier = INaCa_Multiplier;
    extra_data->INaK_Multiplier = INaK_Multiplier;
    extra_data->INab_Multiplier = INab_Multiplier;
    extra_data->ICab_Multiplier = ICab_Multiplier;
    extra_data->IpCa_Multiplier = IpCa_Multiplier;
    extra_data->ICaCl_Multiplier = ICaCl_Multiplier;
    extra_data->IClb_Multiplier = IClb_Multiplier;
    extra_data->Jrel_Multiplier = Jrel_Multiplier;
    extra_data->Jup_Multiplier = Jup_Multiplier;
    
    extra_data->initial_ss_endo = MALLOC_ARRAY_OF_TYPE(real, 43);
    extra_data->initial_ss_epi = MALLOC_ARRAY_OF_TYPE(real, 43);
    extra_data->initial_ss_mid = MALLOC_ARRAY_OF_TYPE(real, 43);

    // Set the initial conditions (celltype = ENDO)
    extra_data->initial_ss_endo[0]  = -8.890585e+01;
    extra_data->initial_ss_endo[1]  = 1.107642e-02;
    extra_data->initial_ss_endo[2]  = 6.504164e-05;
    extra_data->initial_ss_endo[3]  = 1.210818e+01;
    extra_data->initial_ss_endo[4]  = 1.210851e+01;
    extra_data->initial_ss_endo[5]  = 1.426206e+02;
    extra_data->initial_ss_endo[6]  = 1.426205e+02;
    extra_data->initial_ss_endo[7]  = 1.530373e+00;
    extra_data->initial_ss_endo[8]  = 1.528032e+00;
    extra_data->initial_ss_endo[9]  = 7.455488e-05;
    extra_data->initial_ss_endo[10] = 7.814592e-04;
    extra_data->initial_ss_endo[11] = 8.313839e-01;
    extra_data->initial_ss_endo[12] = 8.311938e-01;
    extra_data->initial_ss_endo[13] = 6.752873e-01;
    extra_data->initial_ss_endo[14] = 8.308255e-01;
    extra_data->initial_ss_endo[15] = 1.585610e-04;
    extra_data->initial_ss_endo[16] = 5.294475e-01;
    extra_data->initial_ss_endo[17] = 2.896996e-01;
    extra_data->initial_ss_endo[18] = 9.419166e-04;
    extra_data->initial_ss_endo[19] = 9.996194e-01;
    extra_data->initial_ss_endo[20] = 5.938602e-01;
    extra_data->initial_ss_endo[21] = 4.799180e-04;
    extra_data->initial_ss_endo[22] = 9.996194e-01;
    extra_data->initial_ss_endo[23] = 6.543754e-01;
    extra_data->initial_ss_endo[24] = -2.898677e-33;
    extra_data->initial_ss_endo[25] = 1.000000e+00;
    extra_data->initial_ss_endo[26] = 9.389659e-01;
    extra_data->initial_ss_endo[27] = 1.000000e+00;
    extra_data->initial_ss_endo[28] = 9.999003e-01;
    extra_data->initial_ss_endo[29] = 9.999773e-01;
    extra_data->initial_ss_endo[30] = 1.000000e+00;
    extra_data->initial_ss_endo[31] = 1.000000e+00;
    extra_data->initial_ss_endo[32] = 4.920606e-04;
    extra_data->initial_ss_endo[33] = 8.337021e-04;
    extra_data->initial_ss_endo[34] = 6.962775e-04;
    extra_data->initial_ss_endo[35] = 8.425453e-04;
    extra_data->initial_ss_endo[36] = 9.980807e-01;
    extra_data->initial_ss_endo[37] = 1.289824e-05;
    extra_data->initial_ss_endo[38] = 3.675442e-04;
    extra_data->initial_ss_endo[39] = 2.471690e-01;
    extra_data->initial_ss_endo[40] = 1.742987e-04;
    extra_data->initial_ss_endo[41] = 5.421027e-24;
    extra_data->initial_ss_endo[42] = 6.407933e-23;

    // Set the initial conditions (celltype = EPI)
    extra_data->initial_ss_epi[0]  = -8.917755e+01;
    extra_data->initial_ss_epi[1]  = 1.288116e-02;
    extra_data->initial_ss_epi[2]  = 5.767956e-05;
    extra_data->initial_ss_epi[3]  = 1.284260e+01;
    extra_data->initial_ss_epi[4]  = 1.284291e+01;
    extra_data->initial_ss_epi[5]  = 1.429114e+02;
    extra_data->initial_ss_epi[6]  = 1.429113e+02;
    extra_data->initial_ss_epi[7]  = 1.812268e+00;
    extra_data->initial_ss_epi[8]  = 1.810520e+00;
    extra_data->initial_ss_epi[9]  = 6.631866e-05;
    extra_data->initial_ss_epi[10] = 7.370422e-04;
    extra_data->initial_ss_epi[11] = 8.366816e-01;
    extra_data->initial_ss_epi[12] = 8.366012e-01;
    extra_data->initial_ss_epi[13] = 6.840260e-01;
    extra_data->initial_ss_epi[14] = 8.363958e-01;
    extra_data->initial_ss_epi[15] = 1.505860e-04;
    extra_data->initial_ss_epi[16] = 5.412669e-01;
    extra_data->initial_ss_epi[17] = 3.043382e-01;
    extra_data->initial_ss_epi[18] = 9.248184e-04;
    extra_data->initial_ss_epi[19] = 9.996371e-01;
    extra_data->initial_ss_epi[20] = 9.996342e-01;
    extra_data->initial_ss_epi[21] = 4.712023e-04;
    extra_data->initial_ss_epi[22] = 9.996371e-01;
    extra_data->initial_ss_epi[23] = 9.996366e-01;
    extra_data->initial_ss_epi[24] = 4.333129e-43;
    extra_data->initial_ss_epi[25] = 1.000000e+00;
    extra_data->initial_ss_epi[26] = 9.485160e-01;
    extra_data->initial_ss_epi[27] = 1.000000e+00;
    extra_data->initial_ss_epi[28] = 9.999339e-01;
    extra_data->initial_ss_epi[29] = 9.999822e-01;
    extra_data->initial_ss_epi[30] = 1.000000e+00;
    extra_data->initial_ss_epi[31] = 1.000000e+00;
    extra_data->initial_ss_epi[32] = 3.086885e-04;
    extra_data->initial_ss_epi[33] = 5.303737e-04;
    extra_data->initial_ss_epi[34] = 6.775197e-04;
    extra_data->initial_ss_epi[35] = 8.264829e-04;
    extra_data->initial_ss_epi[36] = 9.982135e-01;
    extra_data->initial_ss_epi[37] = 9.433146e-06;
    extra_data->initial_ss_epi[38] = 2.730221e-04;
    extra_data->initial_ss_epi[39] = 2.308784e-01;
    extra_data->initial_ss_epi[40] = 1.690386e-04;
    extra_data->initial_ss_epi[41] = -1.103286e-23;
    extra_data->initial_ss_epi[42] = -6.177055e-22;

    // Set the initial conditions (celltype = MCELL)
    extra_data->initial_ss_mid[0]  = -8.924177e+01;
    extra_data->initial_ss_mid[1]  = 1.922391e-02;
    extra_data->initial_ss_mid[2]  = 6.585066e-05;
    extra_data->initial_ss_mid[3]  = 1.503347e+01;
    extra_data->initial_ss_mid[4]  = 1.503401e+01;
    extra_data->initial_ss_mid[5]  = 1.434407e+02;
    extra_data->initial_ss_mid[6]  = 1.434406e+02;
    extra_data->initial_ss_mid[7]  = 1.959747e+00;
    extra_data->initial_ss_mid[8]  = 1.963459e+00;
    extra_data->initial_ss_mid[9]  = 8.177438e-05;
    extra_data->initial_ss_mid[10] = 7.269124e-04;
    extra_data->initial_ss_mid[11] = 8.379059e-01;
    extra_data->initial_ss_mid[12] = 8.377164e-01;
    extra_data->initial_ss_mid[13] = 6.860578e-01;
    extra_data->initial_ss_mid[14] = 8.372100e-01;
    extra_data->initial_ss_mid[15] = 1.487602e-04;
    extra_data->initial_ss_mid[16] = 5.350003e-01;
    extra_data->initial_ss_mid[17] = 2.851164e-01;
    extra_data->initial_ss_mid[18] = 9.208259e-04;
    extra_data->initial_ss_mid[19] = 9.996411e-01;
    extra_data->initial_ss_mid[20] = 5.673539e-01;
    extra_data->initial_ss_mid[21] = 4.691672e-04;
    extra_data->initial_ss_mid[22] = 9.996412e-01;
    extra_data->initial_ss_mid[23] = 6.265825e-01;
    extra_data->initial_ss_mid[24] = -4.922960e-40;
    extra_data->initial_ss_mid[25] = 1.000000e+00;
    extra_data->initial_ss_mid[26] = 9.200354e-01;
    extra_data->initial_ss_mid[27] = 1.000000e+00;
    extra_data->initial_ss_mid[28] = 9.997888e-01;
    extra_data->initial_ss_mid[29] = 9.999665e-01;
    extra_data->initial_ss_mid[30] = 1.000000e+00;
    extra_data->initial_ss_mid[31] = 1.000000e+00;
    extra_data->initial_ss_mid[32] = 5.161178e-04;
    extra_data->initial_ss_mid[33] = 1.189422e-03;
    extra_data->initial_ss_mid[34] = 6.917041e-04;
    extra_data->initial_ss_mid[35] = 8.225453e-04;
    extra_data->initial_ss_mid[36] = 9.979358e-01;
    extra_data->initial_ss_mid[37] = 1.835276e-05;
    extra_data->initial_ss_mid[38] = 5.316232e-04;
    extra_data->initial_ss_mid[39] = 2.650323e-01;
    extra_data->initial_ss_mid[40] = 1.678628e-04;
    extra_data->initial_ss_mid[41] = 2.091039e-25;
    extra_data->initial_ss_mid[42] = 2.438403e-23;

    extra_data->transmurality = MALLOC_ARRAY_OF_TYPE(real, num_cells);

    return extra_data;
}

struct extra_data_for_torord * set_common_torord_dyncl_data (struct config *config, uint32_t num_cells) {
    struct extra_data_for_torord *extra_data = MALLOC_ONE_TYPE(struct extra_data_for_torord);

    real INa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INa_Multiplier, config, "INa_Multiplier");

    real ICaL_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaL_Multiplier, config, "ICaL_Multiplier");

    real Ito_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ito_Multiplier, config, "Ito_Multiplier");

    real INaL_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaL_Multiplier, config, "INaL_Multiplier");

    real IKr_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKr_Multiplier, config, "IKr_Multiplier");

    real IKs_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKs_Multiplier, config, "IKs_Multiplier");

    real IK1_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IK1_Multiplier, config, "IK1_Multiplier");

    real IKb_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKb_Multiplier, config, "IKb_Multiplier");

    real INaCa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaCa_Multiplier, config, "INaCa_Multiplier");

    real INaK_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaK_Multiplier, config, "INaK_Multiplier");

    real INab_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INab_Multiplier, config, "INab_Multiplier");

    real ICab_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICab_Multiplier, config, "ICab_Multiplier");

    real IpCa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IpCa_Multiplier, config, "IpCa_Multiplier");

    real ICaCl_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaCl_Multiplier, config, "ICaCl_Multiplier");

    real IClb_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IClb_Multiplier, config, "IClb_Multiplier");

    real Jrel_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Jrel_Multiplier, config, "Jrel_Multiplier");

    real Jup_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Jup_Multiplier, config, "Jup_Multiplier");

    extra_data->INa_Multiplier = INa_Multiplier;
    extra_data->ICaL_Multiplier = ICaL_Multiplier;
    extra_data->Ito_Multiplier = Ito_Multiplier;
    extra_data->INaL_Multiplier = INaL_Multiplier;
    extra_data->IKr_Multiplier = IKr_Multiplier;
    extra_data->IKs_Multiplier = IKs_Multiplier;
    extra_data->IK1_Multiplier = IK1_Multiplier;
    extra_data->IKb_Multiplier = IKb_Multiplier;
    extra_data->INaCa_Multiplier = INaCa_Multiplier;
    extra_data->INaK_Multiplier = INaK_Multiplier;
    extra_data->INab_Multiplier = INab_Multiplier;
    extra_data->ICab_Multiplier = ICab_Multiplier;
    extra_data->IpCa_Multiplier = IpCa_Multiplier;
    extra_data->ICaCl_Multiplier = ICaCl_Multiplier;
    extra_data->IClb_Multiplier = IClb_Multiplier;
    extra_data->Jrel_Multiplier = Jrel_Multiplier;
    extra_data->Jup_Multiplier = Jup_Multiplier;
    
    extra_data->initial_ss_endo = MALLOC_ARRAY_OF_TYPE(real, 45);
    extra_data->initial_ss_epi = MALLOC_ARRAY_OF_TYPE(real, 45);
    extra_data->initial_ss_mid = MALLOC_ARRAY_OF_TYPE(real, 45);

    // Set the initial conditions 200 beats (celltype = ENDO) 
    extra_data->initial_ss_endo[0] = -9.035192e+01;
    extra_data->initial_ss_endo[1] = 1.162900e-02;
    extra_data->initial_ss_endo[2] = 6.500000e-05;
    extra_data->initial_ss_endo[3] = 1.239893e+01;
    extra_data->initial_ss_endo[4] = 1.239926e+01;
    extra_data->initial_ss_endo[5] = 1.482415e+02;
    extra_data->initial_ss_endo[6] = 1.482414e+02;
    extra_data->initial_ss_endo[7] = 1.527292e+00;
    extra_data->initial_ss_endo[8] = 1.524395e+00;
    extra_data->initial_ss_endo[9] = 7.400000e-05;
    extra_data->initial_ss_endo[10] = 5.720000e-04;
    extra_data->initial_ss_endo[11] = 8.579420e-01;
    extra_data->initial_ss_endo[12] = 8.577990e-01;
    extra_data->initial_ss_endo[13] = 7.199660e-01;
    extra_data->initial_ss_endo[14] = 8.575760e-01;
    extra_data->initial_ss_endo[15] = 1.200000e-04;
    extra_data->initial_ss_endo[16] = 5.748970e-01;
    extra_data->initial_ss_endo[17] = 3.250180e-01;
    extra_data->initial_ss_endo[18] = 8.540000e-04;
    extra_data->initial_ss_endo[19] = 9.997050e-01;
    extra_data->initial_ss_endo[20] = 5.959350e-01;
    extra_data->initial_ss_endo[21] = 4.350000e-04;
    extra_data->initial_ss_endo[22] = 9.997050e-01;
    extra_data->initial_ss_endo[23] = 6.589890e-01;
    extra_data->initial_ss_endo[24] = 0.000000e+00;
    extra_data->initial_ss_endo[25] = 1.000000e+00;
    extra_data->initial_ss_endo[26] = 9.343710e-01;
    extra_data->initial_ss_endo[27] = 1.000000e+00;
    extra_data->initial_ss_endo[28] = 9.998810e-01;
    extra_data->initial_ss_endo[29] = 9.999820e-01;
    extra_data->initial_ss_endo[30] = 1.000000e+00;
    extra_data->initial_ss_endo[31] = 1.000000e+00;
    extra_data->initial_ss_endo[32] = 4.830000e-04;
    extra_data->initial_ss_endo[33] = 8.180000e-04;
    extra_data->initial_ss_endo[34] = 9.983340e-01;
    extra_data->initial_ss_endo[35] = 7.600000e-04;
    extra_data->initial_ss_endo[36] = 6.260000e-04;
    extra_data->initial_ss_endo[37] = 9.000000e-06;
    extra_data->initial_ss_endo[38] = 2.720000e-04;
    extra_data->initial_ss_endo[39] = 2.568150e-01;
    extra_data->initial_ss_endo[40] = 1.480000e-04;
    extra_data->initial_ss_endo[41] = 0.000000e+00;
    extra_data->initial_ss_endo[42] = 0.000000e+00;
    extra_data->initial_ss_endo[43] = 2.978204e+01;
    extra_data->initial_ss_endo[44] = 2.978201e+01;

    // Set the initial conditions 200 beats (celltype = EPI)
    extra_data->initial_ss_epi[0] = -9.038446e+01;
    extra_data->initial_ss_epi[1] = 1.321400e-02;
    extra_data->initial_ss_epi[2] = 5.700000e-05;
    extra_data->initial_ss_epi[3] = 1.296837e+01;
    extra_data->initial_ss_epi[4] = 1.296867e+01;
    extra_data->initial_ss_epi[5] = 1.477895e+02;
    extra_data->initial_ss_epi[6] = 1.477895e+02;
    extra_data->initial_ss_epi[7] = 1.801231e+00;
    extra_data->initial_ss_epi[8] = 1.798984e+00;
    extra_data->initial_ss_epi[9] = 6.600000e-05;
    extra_data->initial_ss_epi[10] = 5.680000e-04;
    extra_data->initial_ss_epi[11] = 8.585000e-01;
    extra_data->initial_ss_epi[12] = 8.584370e-01;
    extra_data->initial_ss_epi[13] = 7.209280e-01;
    extra_data->initial_ss_epi[14] = 8.583060e-01;
    extra_data->initial_ss_epi[15] = 1.200000e-04;
    extra_data->initial_ss_epi[16] = 5.793440e-01;
    extra_data->initial_ss_epi[17] = 3.348040e-01;
    extra_data->initial_ss_epi[18] = 8.530000e-04;
    extra_data->initial_ss_epi[19] = 9.997060e-01;
    extra_data->initial_ss_epi[20] = 9.997050e-01;
    extra_data->initial_ss_epi[21] = 4.340000e-04;
    extra_data->initial_ss_epi[22] = 9.997060e-01;
    extra_data->initial_ss_epi[23] = 9.997060e-01;
    extra_data->initial_ss_epi[24] = 0.000000e+00;
    extra_data->initial_ss_epi[25] = 1.000000e+00;
    extra_data->initial_ss_epi[26] = 9.436160e-01;
    extra_data->initial_ss_epi[27] = 1.000000e+00;
    extra_data->initial_ss_epi[28] = 9.999210e-01;
    extra_data->initial_ss_epi[29] = 9.999850e-01;
    extra_data->initial_ss_epi[30] = 1.000000e+00;
    extra_data->initial_ss_epi[31] = 1.000000e+00;
    extra_data->initial_ss_epi[32] = 3.020000e-04;
    extra_data->initial_ss_epi[33] = 5.160000e-04;
    extra_data->initial_ss_epi[34] = 9.984000e-01;
    extra_data->initial_ss_epi[35] = 7.590000e-04;
    extra_data->initial_ss_epi[36] = 6.200000e-04;
    extra_data->initial_ss_epi[37] = 7.000000e-06;
    extra_data->initial_ss_epi[38] = 2.140000e-04;
    extra_data->initial_ss_epi[39] = 2.426290e-01;
    extra_data->initial_ss_epi[40] = 1.480000e-04;
    extra_data->initial_ss_epi[41] = 0.000000e+00;
    extra_data->initial_ss_epi[42] = 0.000000e+00;
    extra_data->initial_ss_epi[43] = 3.005381e+01;
    extra_data->initial_ss_epi[44] = 3.005378e+01;

    // Set the initial conditions 200 beats (celltype = MCELL)
    extra_data->initial_ss_mid[0] = -9.032701e+01;
    extra_data->initial_ss_mid[1] = 1.889700e-02;
    extra_data->initial_ss_mid[2] = 6.400000e-05;
    extra_data->initial_ss_mid[3] = 1.334108e+01;
    extra_data->initial_ss_mid[4] = 1.334155e+01;
    extra_data->initial_ss_mid[5] = 1.474754e+02;
    extra_data->initial_ss_mid[6] = 1.474753e+02;
    extra_data->initial_ss_mid[7] = 1.814771e+00;
    extra_data->initial_ss_mid[8] = 1.815495e+00;
    extra_data->initial_ss_mid[9] = 7.800000e-05;
    extra_data->initial_ss_mid[10] = 5.750000e-04;
    extra_data->initial_ss_mid[11] = 8.575170e-01;
    extra_data->initial_ss_mid[12] = 8.573630e-01;
    extra_data->initial_ss_mid[13] = 7.192340e-01;
    extra_data->initial_ss_mid[14] = 8.569750e-01;
    extra_data->initial_ss_mid[15] = 1.210000e-04;
    extra_data->initial_ss_mid[16] = 5.666070e-01;
    extra_data->initial_ss_mid[17] = 3.061970e-01;
    extra_data->initial_ss_mid[18] = 8.560000e-04;
    extra_data->initial_ss_mid[19] = 9.997030e-01;
    extra_data->initial_ss_mid[20] = 5.577600e-01;
    extra_data->initial_ss_mid[21] = 4.360000e-04;
    extra_data->initial_ss_mid[22] = 9.997030e-01;
    extra_data->initial_ss_mid[23] = 6.157570e-01;
    extra_data->initial_ss_mid[24] = 0.000000e+00;
    extra_data->initial_ss_mid[25] = 1.000000e+00;
    extra_data->initial_ss_mid[26] = 9.038390e-01;
    extra_data->initial_ss_mid[27] = 1.000000e+00;
    extra_data->initial_ss_mid[28] = 9.996760e-01;
    extra_data->initial_ss_mid[29] = 9.999600e-01;
    extra_data->initial_ss_mid[30] = 1.000000e+00;
    extra_data->initial_ss_mid[31] = 1.000000e+00;
    extra_data->initial_ss_mid[32] = 4.620000e-04;
    extra_data->initial_ss_mid[33] = 9.760000e-04;
    extra_data->initial_ss_mid[34] = 9.980740e-01;
    extra_data->initial_ss_mid[35] = 7.610000e-04;
    extra_data->initial_ss_mid[36] = 6.430000e-04;
    extra_data->initial_ss_mid[37] = 1.600000e-05;
    extra_data->initial_ss_mid[38] = 5.050000e-04;
    extra_data->initial_ss_mid[39] = 3.021660e-01;
    extra_data->initial_ss_mid[40] = 1.490000e-04;
    extra_data->initial_ss_mid[41] = 0.000000e+00;
    extra_data->initial_ss_mid[42] = 0.000000e+00;
    extra_data->initial_ss_mid[43] = 3.026206e+01;
    extra_data->initial_ss_mid[44] = 3.026203e+01;

    extra_data->transmurality = MALLOC_ARRAY_OF_TYPE(real, num_cells);

    return extra_data;
}

struct extra_data_for_trovato * set_common_trovato_data (struct config *config, uint32_t num_cells) {
    struct extra_data_for_trovato *extra_data = MALLOC_ONE_TYPE(struct extra_data_for_trovato);

    real GNa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GNa_Multiplier, config, "GNa_Multiplier");
    real GNaL_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GNaL_Multiplier, config, "GNaL_Multiplier");
    real GCaT_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GCaT_Multiplier, config, "GCaT_Multiplier");
    real Gto_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Gto_Multiplier, config, "Gto_Multiplier");
    real Gsus_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Gsus_Multiplier, config, "Gsus_Multiplier");
    real Gkr_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Gkr_Multiplier, config, "Gkr_Multiplier");
    real Gks_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Gks_Multiplier, config, "Gks_Multiplier");
    real GfNa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GfNa_Multiplier, config, "GfNa_Multiplier");
    real GfK_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GfK_Multiplier, config, "GfK_Multiplier");
    real GK1_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GK1_Multiplier, config, "GK1_Multiplier");
    real GNCX_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GNCX_Multiplier, config, "GNCX_Multiplier");
    real GNaK_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GNaK_Multiplier, config, "GNaK_Multiplier");
    real INa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INa_Multiplier, config, "INa_Multiplier");
    real ICaL_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaL_Multiplier, config, "ICaL_Multiplier");
    real ICaNa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaNa_Multiplier, config, "ICaNa_Multiplier");
    real ICaK_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaK_Multiplier, config, "ICaK_Multiplier");
    real Ito_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ito_Multiplier, config, "Ito_Multiplier");
    real INaL_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaL_Multiplier, config, "INaL_Multiplier");
    real IKr_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKr_Multiplier, config, "IKr_Multiplier");
    real IKs_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKs_Multiplier, config, "IKs_Multiplier");
    real IK1_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IK1_Multiplier, config, "IK1_Multiplier");
    real INaCa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaCa_Multiplier, config, "INaCa_Multiplier");
    real INaK_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaK_Multiplier, config, "INaK_Multiplier");
    real INab_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INab_Multiplier, config, "INab_Multiplier");
    real ICab_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICab_Multiplier, config, "ICab_Multiplier");
    real ICaT_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaT_Multiplier, config, "ICaT_Multiplier");
    real Isus_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Isus_Multiplier, config, "Isus_Multiplier");
    real If_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, If_Multiplier, config, "If_Multiplier");
    real IpCa_Multiplier = 1.0;    
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IpCa_Multiplier, config, "IpCa_Multiplier");

    extra_data->GNa_Multiplier  = GNa_Multiplier;
    extra_data->GNaL_Multiplier = GNaL_Multiplier;
    extra_data->GCaT_Multiplier = GCaT_Multiplier;
    extra_data->Gto_Multiplier  = Gto_Multiplier;
    extra_data->Gsus_Multiplier = Gsus_Multiplier;
    extra_data->Gkr_Multiplier  = Gkr_Multiplier;
    extra_data->Gks_Multiplier  = Gks_Multiplier;
    extra_data->GfNa_Multiplier = GfNa_Multiplier;
    extra_data->GfK_Multiplier  = GfK_Multiplier;
    extra_data->GK1_Multiplier  = GK1_Multiplier;
    extra_data->GNCX_Multiplier = GNCX_Multiplier;
    extra_data->GNaK_Multiplier = GNaK_Multiplier;
    extra_data->INa_Multiplier  = INa_Multiplier; 
    extra_data->ICaL_Multiplier = ICaL_Multiplier;
    extra_data->ICaNa_Multiplier= ICaNa_Multiplier;
    extra_data->ICaK_Multiplier = ICaK_Multiplier;
    extra_data->Ito_Multiplier  = Ito_Multiplier;
    extra_data->INaL_Multiplier = INaL_Multiplier;
    extra_data->IKr_Multiplier  = IKr_Multiplier; 
    extra_data->IKs_Multiplier  = IKs_Multiplier; 
    extra_data->IK1_Multiplier  = IK1_Multiplier; 
    extra_data->INaCa_Multiplier= INaCa_Multiplier;
    extra_data->INaK_Multiplier = INaK_Multiplier;  
    extra_data->INab_Multiplier = INab_Multiplier;  
    extra_data->ICab_Multiplier = ICab_Multiplier;
    extra_data->ICaT_Multiplier = ICaT_Multiplier;
    extra_data->Isus_Multiplier = Isus_Multiplier;
    extra_data->If_Multiplier   = If_Multiplier;
    extra_data->IpCa_Multiplier = IpCa_Multiplier;

    return extra_data;
}