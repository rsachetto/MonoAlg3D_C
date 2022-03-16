#include "../libraries_common/common_data_structures.h"
#include <stdio.h>

struct pseudo_bidomain_persistent_data {
     real_cpu **distances;
     real_cpu *main_diagonal;
     real_cpu *beta_im;
     struct point_3d *leads;
     FILE *output_file;
     real_cpu scale_factor;
     int n_leads;
};   

#define PSEUDO_BIDOMAIN_DATA ((struct pseudo_bidomain_persistent_data *)config->persistent_data)

