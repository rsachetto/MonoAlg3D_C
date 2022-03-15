//
// Created by sachetto on 15/03/22.
//

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"
#include "../utils/utils.h"
#include "../config/ecg_config.h"
#include <time.h>
#include <unistd.h>

INIT_CALC_ECG(init_pseudo_bidomain) {
    log_info("INIT PSEUDO ECG\n");
}

CALC_ECG(pseudo_bidomain) {
    log_info("CALC PSEUDO ECG\n");
}

INIT_CALC_ECG(end_pseudo_bidomain) {
    log_info("END PSEUDO ECG\n");
}

