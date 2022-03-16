//
// Created by sachetto on 15/03/22.
//

#include "ecg.h"
#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../config/ecg_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"
#include "../utils/utils.h"
#include <time.h>
#include <unistd.h>
#include <float.h>

#define EUCLIDIAN_DISTANCE(p, q) sqrt(pow(p.x - q.x, 2.0) + pow(p.y - q.y, 2.0) + pow(p.z - q.z, 2.0))

static void get_leads(struct config *config, struct pseudo_bidomain_persistent_data *data) {

    char lead_name[1024];
    int count = 1;

    real_cpu leads_coords[3] = {FLT_MAX,FLT_MAX, FLT_MAX};

    while(true) {
        sprintf(lead_name, "lead%d", count);

        GET_PARAMETER_VECTOR3_VALUE_OR_USE_DEFAULT(leads_coords, config, lead_name);

        if(count == 1 && leads_coords[0] == FLT_MAX) {
            log_error_and_exit("No leads defined on [calc_ecg]!\n");
        }

        if(leads_coords[0] == FLT_MAX) {
            break;
        }

        struct point_3d coord = POINT3D(leads_coords[0], leads_coords[1], leads_coords[2]);
        arrpush(PSEUDO_BIDOMAIN_DATA->leads, coord);

        leads_coords[0] = leads_coords[1] = leads_coords[2] = FLT_MAX;

        count++;
    }
}

INIT_CALC_ECG(init_pseudo_bidomain) {

    config->persistent_data = CALLOC_ONE_TYPE(struct pseudo_bidomain_persistent_data);

    char *filename = strdup("./ecg.txt");
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(filename, config, "filename");
   
    PSEUDO_BIDOMAIN_DATA->output_file = fopen(filename, "w");

    if(PSEUDO_BIDOMAIN_DATA->output_file == NULL) {
        log_error_and_exit("init_pseudo_bidomain - Unable to open file %s!\n", filename);
    }

    real_cpu sigma_b = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_b, config, "sigma_b");

    if(sigma_b == 0.0) {
        log_error_and_exit("init_pseudo_bidomain - sigma_b can't be 0!\n");
    }

    PSEUDO_BIDOMAIN_DATA->scale_factor = 1.0/(4.0*M_PI*sigma_b);

    get_leads(config, PSEUDO_BIDOMAIN_DATA);
    PSEUDO_BIDOMAIN_DATA->n_leads = arrlen(PSEUDO_BIDOMAIN_DATA->leads);

    int n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    PSEUDO_BIDOMAIN_DATA->distances = MALLOC_ARRAY_OF_TYPE(real_cpu *, PSEUDO_BIDOMAIN_DATA->n_leads);
    for(int i = 0; i < PSEUDO_BIDOMAIN_DATA->n_leads; i++) {
        PSEUDO_BIDOMAIN_DATA->distances[i] = MALLOC_ARRAY_OF_TYPE(real_cpu, n_active);
    }

    PSEUDO_BIDOMAIN_DATA->beta_im = MALLOC_ARRAY_OF_TYPE(real_cpu, n_active);

    PSEUDO_BIDOMAIN_DATA->main_diagonal = MALLOC_ARRAY_OF_TYPE(real_cpu, n_active);

    // calc the distances from each volume to each electrode (r)
    for(int i = 0; i < PSEUDO_BIDOMAIN_DATA->n_leads; i++) {

        struct point_3d lead = PSEUDO_BIDOMAIN_DATA->leads[i];

        OMP(parallel for)
        for(int j = 0; j < n_active; j++) {
            struct point_3d center = ac[j]->center;
            PSEUDO_BIDOMAIN_DATA->distances[i][j] = EUCLIDIAN_DISTANCE(lead, center);
        }
    }

    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    // calc the main diagonal for the ECG calculation (matrix main diag - ALPHA)
    OMP(parallel for)
    for(int i = 0; i < n_active; i++) {
        real_cpu alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        PSEUDO_BIDOMAIN_DATA->main_diagonal[i] = ac[i]->elements[0].value - alpha;
    }

    free(filename);
}

CALC_ECG(pseudo_bidomain) {
    // log_info("CALC PSEUDO ECG\n");
    // use the equation described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378475/#FD7
    // for each electrode, integrate ECG Matrix x Vm (maybe we can save Vm to avoid extra GPU copy)
    
    int n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    OMP(parallel for)
    for(int i = 0; i < n_active; i++) {

        struct element *cell_elements = ac[i]->elements;
        size_t max_el = arrlen(cell_elements);

        PSEUDO_BIDOMAIN_DATA->beta_im[i] = PSEUDO_BIDOMAIN_DATA->main_diagonal[i] * cell_elements[0].cell->v;

        for(size_t el = 1; el < max_el; el++) {
            PSEUDO_BIDOMAIN_DATA->beta_im[i] += cell_elements[el].value * cell_elements[el].cell->v;
        }
    }

    fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "%lf ", time_info->current_t);

    for(int i = 0; i < PSEUDO_BIDOMAIN_DATA->n_leads; i++) {
        real_cpu local_sum = 0;
        
        OMP(parallel for reduction(+:local_sum))
        for(int j = 0; j < n_active; j++) {
            struct point_3d d = ac[j]->discretization;
            local_sum += (PSEUDO_BIDOMAIN_DATA->beta_im[j] / PSEUDO_BIDOMAIN_DATA->distances[i][j]) * (d.x*d.y*d.z);
        }

        fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "%lf ", PSEUDO_BIDOMAIN_DATA->scale_factor*local_sum);
    }

    fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "\n");
}

INIT_CALC_ECG(end_pseudo_bidomain) {
    // Free distances and maybe the saved Vm
    fclose(PSEUDO_BIDOMAIN_DATA->output_file);
}
