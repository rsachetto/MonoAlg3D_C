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

    //calc the distances from each volume to each electrode (r)
    //calc the main diagonal for the ECG calculation (matrix main diag - ALPHA)

}

CALC_ECG(pseudo_bidomain) {
    log_info("CALC PSEUDO ECG\n");
    //use the equation described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378475/#FD7
    //for each electrode, integrate ECG Matrix x Vm (maybe we can save Vm to avoid extra GPU copy)
}

INIT_CALC_ECG(end_pseudo_bidomain) {
    log_info("END PSEUDO ECG\n");
    //Free distances and maybe the saved Vm
}

