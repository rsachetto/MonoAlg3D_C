//
// Created by sachetto on 15/03/22.
//

#include "ecg.h"
#include "../3dparty/stb_ds.h"
#include "../config/ecg_config.h"
#include "../utils/file_utils.h"
#include "../utils/utils.h"
#include <float.h>

#ifdef COMPILE_CUDA

#include "../gpu_utils/gpu_utils.h"

// Precision to be used for the calculations on the GPU
#ifdef CELL_MODEL_REAL_DOUBLE
#pragma message("calc_ecg, using double precision on the GPU")
#define CUBLAS_SIZE CUDA_R_64F
#else
#pragma message("calc_ecg, using single precision on the GPU")
#define CUBLAS_SIZE CUDA_R_32F
#endif

#if CUSPARSE_VER_MAJOR >= 12
#define CUSPARSE_ALG CUSPARSE_SPMV_ALG_DEFAULT
#elif CUSPARSE_VER_MAJOR == 11
#define CUSPARSE_ALG CUSPARSE_MV_ALG_DEFAULT
#endif

#endif // COMPILE_CUDA

#define EUCLIDIAN_DISTANCE(p, q) sqrt(pow(p.x - q.x, 2.0) + pow(p.y - q.y, 2.0) + pow(p.z - q.z, 2.0))

void fill_element(uint32_t position, real_cpu dx, real_cpu dy, real_cpu dz, struct element *cell_elements) {

    real_cpu multiplier;
    multiplier = ((dx * dy) / dz);
    cell_elements[position].value_ecg = -multiplier;
    cell_elements[0].value_ecg += (multiplier);
}

static void fill_discretization_matrix_elements(struct cell_node *grid_cell, void *neighbour_grid_cell, enum transition_direction direction) {

    (void) direction;

    bool has_found;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
* cell which is a cell node. */
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    enum cell_type neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if(neighbour_grid_cell_level > grid_cell->cell_data.level) {
        if(neighbour_grid_cell_type == TRANSITION_NODE) {
            has_found = false;
            while(!has_found) {
                if(neighbour_grid_cell_type == TRANSITION_NODE) {
                    white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                    if(white_neighbor_cell->single_connector == NULL) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    } else {
        if(neighbour_grid_cell_level <= grid_cell->cell_data.level && (neighbour_grid_cell_type == TRANSITION_NODE)) {
            has_found = false;
            while(!has_found) {
                if(neighbour_grid_cell_type == TRANSITION_NODE) {
                    white_neighbor_cell = (struct transition_node *)(neighbour_grid_cell);
                    if(white_neighbor_cell->single_connector == 0) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    }

    // We care only with the interior points
    if(neighbour_grid_cell_type == CELL_NODE) {

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if(black_neighbor_cell->active) {

            uint32_t position;
            real_cpu dx, dy, dz;

            if(black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) {
                dx = black_neighbor_cell->discretization.x;
                dy = black_neighbor_cell->discretization.y;
                dz = black_neighbor_cell->discretization.z;
            } else {
                dx = grid_cell->discretization.x;
                dy = grid_cell->discretization.y;
                dz = grid_cell->discretization.z;
            }

            lock_cell_node(grid_cell);

            struct element *cell_elements = grid_cell->elements;
            position = black_neighbor_cell->grid_position;

            size_t max_elements = arrlen(cell_elements);

        bool insert = false;
        size_t p;

        for(size_t i = 0; i < max_elements; i++) {
            if(cell_elements[i].column == position) {
            p = i;
            insert = true;
                break;
            }
        }

        if(insert) {
            fill_element(p, dx, dy, dz, cell_elements);
        }
            unlock_cell_node(grid_cell);
        }
    }
}

void assembly_divergent(struct grid *the_grid) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    for(uint32_t i = 0; i < num_active_cells; i++) {

        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }

}

static void get_leads(struct config *config) {

    char lead_name[1024];
    int count = 1;

    real_cpu leads_coords[3] = {FLT_MAX, FLT_MAX, FLT_MAX};

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

INIT_CALC_ECG(init_pseudo_bidomain_cpu) {
    config->persistent_data = CALLOC_ONE_TYPE(struct pseudo_bidomain_persistent_data);

    char *filename = strdup("./ecg.txt");
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(filename, config, "filename");

    char *dir = get_dir_from_path(filename);
    create_dir(dir);
    free(dir);

    PSEUDO_BIDOMAIN_DATA->output_file = fopen(filename, "w");

    if(PSEUDO_BIDOMAIN_DATA->output_file == NULL) {
        log_error_and_exit("init_pseudo_bidomain - Unable to open file %s!\n", filename);
    }

    real_cpu sigma_b = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_b, config, "sigma_b");

    if(sigma_b == 0.0) {
        log_error_and_exit("init_pseudo_bidomain - sigma_b can't be 0!\n");
    }

    PSEUDO_BIDOMAIN_DATA->scale_factor = 1.0 / (4.0 * M_PI * sigma_b);

    get_leads(config);
    PSEUDO_BIDOMAIN_DATA->n_leads = arrlen(PSEUDO_BIDOMAIN_DATA->leads);

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    PSEUDO_BIDOMAIN_DATA->distances = MALLOC_ARRAY_OF_TYPE(real, PSEUDO_BIDOMAIN_DATA->n_leads * n_active);

    PSEUDO_BIDOMAIN_DATA->beta_im = MALLOC_ARRAY_OF_TYPE(real, n_active);

    // calc the distances from each volume to each electrode (r)
    for(uint32_t i = 0; i < PSEUDO_BIDOMAIN_DATA->n_leads; i++) {

        struct point_3d lead = PSEUDO_BIDOMAIN_DATA->leads[i];

        OMP(parallel for)
        for(int j = 0; j < n_active; j++) {
            uint32_t index = i * n_active + j;
            struct point_3d center = ac[j]->center;
            PSEUDO_BIDOMAIN_DATA->distances[index] = EUCLIDIAN_DISTANCE(lead, center);
        }
    }

    assembly_divergent(the_grid);

    free(filename);
}

CALC_ECG(pseudo_bidomain_cpu) {
    // use the equation described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3378475/#FD7
    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    OMP(parallel for)
    for(uint32_t i = 0; i < n_active; i++) {
        struct element *cell_elements = ac[i]->elements;
        size_t max_el = arrlen(cell_elements);

        struct point_3d d = ac[i]->discretization;
        real_cpu volume = d.x * d.y * d.z;

        PSEUDO_BIDOMAIN_DATA->beta_im[i] = 0;

        for(size_t el = 0; el < max_el; el++) {
            PSEUDO_BIDOMAIN_DATA->beta_im[i] += cell_elements[el].value_ecg * cell_elements[el].cell->v / volume;
        }
    }

    fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "%lf ", time_info->current_t);

    for(uint32_t i = 0; i < PSEUDO_BIDOMAIN_DATA->n_leads; i++) {
        real_cpu local_sum = 0.0;

        OMP(parallel for reduction(+:local_sum))
        for(uint32_t j = 0; j < n_active; j++) {
            struct point_3d d = ac[j]->discretization;
            real_cpu volume = d.x * d.y * d.z;

            uint32_t index = i * n_active + j;
            local_sum += ((PSEUDO_BIDOMAIN_DATA->beta_im[j] / PSEUDO_BIDOMAIN_DATA->distances[index])) * volume;
        }

        fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "%lf ", -PSEUDO_BIDOMAIN_DATA->scale_factor * local_sum);
    }

    fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "\n");
}

END_CALC_ECG(end_pseudo_bidomain_cpu) {
    fclose(PSEUDO_BIDOMAIN_DATA->output_file);
    free(PSEUDO_BIDOMAIN_DATA->beta_im);
    arrfree(PSEUDO_BIDOMAIN_DATA->leads);
}

#ifdef COMPILE_CUDA

INIT_CALC_ECG(init_pseudo_bidomain_gpu) {

    init_pseudo_bidomain_cpu(config, NULL, the_grid);

    if(!the_ode_solver->gpu) {
        log_warn("The current implementation of pseudo_bidomain_gpu only works when the odes are also being solved using the GPU! Falling back to CPU version\n");
        shput(config->config_data, "use_gpu", "false");
        return;
    }

    //This is allocated when using the CPU code, but we do not need it in the gpu version
    free(PSEUDO_BIDOMAIN_DATA->beta_im);

    check_cublas_error(cusparseCreate(&(PSEUDO_BIDOMAIN_DATA->cusparseHandle)));
    check_cublas_error(cublasCreate(&(PSEUDO_BIDOMAIN_DATA->cublasHandle)));

    int_array I = NULL, J = NULL;
    f32_array val = NULL;

    grid_to_csr_for_ecg(the_grid, &val, &I, &J, false, true);

    int nz = arrlen(val);
    uint32_t N = the_grid->num_active_cells;

    // We need to do that because val is originally a float value. Maybe we should make it real.
    real *new_val = NULL;
    arrsetlen(new_val, nz);
    for(int i = 0; i < nz; i++) {
        new_val[i] = (real)val[i];
    }

    PSEUDO_BIDOMAIN_DATA->nz = nz;

    check_cuda_error(cudaMalloc((void **)&(PSEUDO_BIDOMAIN_DATA->d_col), nz * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(PSEUDO_BIDOMAIN_DATA->d_row), (N + 1) * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(PSEUDO_BIDOMAIN_DATA->d_val), nz * sizeof(real)));
    check_cuda_error(cudaMalloc((void **)&(PSEUDO_BIDOMAIN_DATA->beta_im), N * sizeof(real)));
    check_cuda_error(cudaMalloc((void **)&(PSEUDO_BIDOMAIN_DATA->d_distances), PSEUDO_BIDOMAIN_DATA->n_leads * N * sizeof(real)));
    check_cuda_error(cudaMalloc((void **)&(PSEUDO_BIDOMAIN_DATA->d_volumes), N * sizeof(real)));
    check_cuda_error(cudaMalloc((void **)&(PSEUDO_BIDOMAIN_DATA->tmp_data), N * sizeof(real)));

#if CUBLAS_VER_MAJOR >= 11
    check_cuda_error(cusparseCreateCsr(&(PSEUDO_BIDOMAIN_DATA->matA), N, N, nz, PSEUDO_BIDOMAIN_DATA->d_row, PSEUDO_BIDOMAIN_DATA->d_col,
                                       PSEUDO_BIDOMAIN_DATA->d_val, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUBLAS_SIZE));
    check_cuda_error(cusparseCreateDnVec(&(PSEUDO_BIDOMAIN_DATA->vec_beta_im), N, PSEUDO_BIDOMAIN_DATA->beta_im, CUBLAS_SIZE));
    check_cuda_error(cusparseCreateDnVec(&(PSEUDO_BIDOMAIN_DATA->vec_vm), N, the_ode_solver->sv, CUBLAS_SIZE));
#else
    check_cuda_error((cudaError_t)cusparseCreateMatDescr(&(PSEUDO_BIDOMAIN_DATA->descr)));
    cusparseSetMatType(PSEUDO_BIDOMAIN_DATA->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(PSEUDO_BIDOMAIN_DATA->descr, CUSPARSE_INDEX_BASE_ZERO);
    PSEUDO_BIDOMAIN_DATA->local_sv = the_ode_solver->sv;
#endif

    check_cuda_error(cudaMemcpy(PSEUDO_BIDOMAIN_DATA->d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice));        // JA
    check_cuda_error(cudaMemcpy(PSEUDO_BIDOMAIN_DATA->d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice));   // IA
    check_cuda_error(cudaMemcpy(PSEUDO_BIDOMAIN_DATA->d_val, new_val, nz * sizeof(real), cudaMemcpyHostToDevice)); // A
    check_cuda_error(cudaMemcpy(PSEUDO_BIDOMAIN_DATA->d_distances, PSEUDO_BIDOMAIN_DATA->distances, PSEUDO_BIDOMAIN_DATA->n_leads * N * sizeof(real),
                                cudaMemcpyHostToDevice));

    PSEUDO_BIDOMAIN_DATA->volumes = MALLOC_ARRAY_OF_TYPE(real, N);

    struct cell_node **ac = the_grid->active_cells;
    OMP(parallel for)
    for(int i = 0; i < N; i += 1) {
        struct point_3d d = ac[i]->discretization;
        PSEUDO_BIDOMAIN_DATA->volumes[i] = d.x * d.y * d.z;
    }

    check_cuda_error(cudaMemcpy(PSEUDO_BIDOMAIN_DATA->d_volumes, PSEUDO_BIDOMAIN_DATA->volumes, N * sizeof(real), cudaMemcpyHostToDevice));

#if CUSPARSE_VER_MAJOR >= 11
    real alpha = 1.0;
    real beta = 0.0;

    check_cuda_error(cusparseSpMV_bufferSize(PSEUDO_BIDOMAIN_DATA->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, PSEUDO_BIDOMAIN_DATA->matA,
                                             PSEUDO_BIDOMAIN_DATA->vec_vm, &beta, PSEUDO_BIDOMAIN_DATA->vec_beta_im, CUBLAS_SIZE, CUSPARSE_ALG,
                                             &(PSEUDO_BIDOMAIN_DATA->bufferSize)));

    check_cuda_error(cudaMalloc(&(PSEUDO_BIDOMAIN_DATA->buffer), PSEUDO_BIDOMAIN_DATA->bufferSize));
#endif

    arrfree(I);
    arrfree(J);
    arrfree(val);
    arrfree(new_val);

    free(PSEUDO_BIDOMAIN_DATA->volumes);
    free(PSEUDO_BIDOMAIN_DATA->distances);
}

CALC_ECG(pseudo_bidomain_gpu) {

    uint32_t n_active = the_grid->num_active_cells;

    // VM is correct
    real alpha = 1.0;
    real beta = 0.0;
#if CUSPARSE_VER_MAJOR >= 11
    check_cublas_error(cusparseSpMV(PSEUDO_BIDOMAIN_DATA->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, PSEUDO_BIDOMAIN_DATA->matA,
                                    PSEUDO_BIDOMAIN_DATA->vec_vm, &beta, PSEUDO_BIDOMAIN_DATA->vec_beta_im, CUBLAS_SIZE, CUSPARSE_ALG,
                                    PSEUDO_BIDOMAIN_DATA->buffer));
#else

#ifdef CELL_MODEL_REAL_DOUBLE
    cusparseDcsrmv(PSEUDO_BIDOMAIN_DATA->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, n_active, n_active, PSEUDO_BIDOMAIN_DATA->nz, &alpha,
                   PSEUDO_BIDOMAIN_DATA->descr, PSEUDO_BIDOMAIN_DATA->d_val, PSEUDO_BIDOMAIN_DATA->d_row, PSEUDO_BIDOMAIN_DATA->d_col,
                   PSEUDO_BIDOMAIN_DATA->local_sv, &beta, PSEUDO_BIDOMAIN_DATA->beta_im);
#else
    cusparseScsrmv(PSEUDO_BIDOMAIN_DATA->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, n_active, n_active, PSEUDO_BIDOMAIN_DATA->nz, &alpha,
                   PSEUDO_BIDOMAIN_DATA->descr, PSEUDO_BIDOMAIN_DATA->d_val, PSEUDO_BIDOMAIN_DATA->d_row, PSEUDO_BIDOMAIN_DATA->d_col,
                   PSEUDO_BIDOMAIN_DATA->local_sv, &beta, PSEUDO_BIDOMAIN_DATA->beta_im);
#endif

#endif

    fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "%lf ", time_info->current_t);

    gpu_vec_div_vec(PSEUDO_BIDOMAIN_DATA->beta_im, PSEUDO_BIDOMAIN_DATA->d_volumes, PSEUDO_BIDOMAIN_DATA->beta_im, n_active);

    for(int i = 0; i < PSEUDO_BIDOMAIN_DATA->n_leads; i++) {

        // beta_im / distance
        gpu_vec_div_vec(PSEUDO_BIDOMAIN_DATA->beta_im, PSEUDO_BIDOMAIN_DATA->d_distances + n_active * i, PSEUDO_BIDOMAIN_DATA->tmp_data, n_active);
        real local_sum;

#ifdef CELL_MODEL_REAL_DOUBLE
        check_cublas_error(
            cublasDdot(PSEUDO_BIDOMAIN_DATA->cublasHandle, n_active, PSEUDO_BIDOMAIN_DATA->tmp_data, 1, PSEUDO_BIDOMAIN_DATA->d_volumes, 1, &local_sum));
#else
        check_cublas_error(
            cublasSdot(PSEUDO_BIDOMAIN_DATA->cublasHandle, n_active, PSEUDO_BIDOMAIN_DATA->tmp_data, 1, PSEUDO_BIDOMAIN_DATA->d_volumes, 1, &local_sum));
#endif

        fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "%lf ", -PSEUDO_BIDOMAIN_DATA->scale_factor * local_sum);
    }

    fprintf(PSEUDO_BIDOMAIN_DATA->output_file, "\n");

}

END_CALC_ECG(end_pseudo_bidomain_gpu) {

    struct pseudo_bidomain_persistent_data *persistent_data = (struct pseudo_bidomain_persistent_data *)config->persistent_data;

    if(!persistent_data)
        return;

    check_cuda_error((cudaError_t)cusparseDestroy(persistent_data->cusparseHandle));
    check_cuda_error((cudaError_t)cublasDestroy(persistent_data->cublasHandle));

#if CUBLAS_VER_MAJOR >= 11
    if(persistent_data->matA) {
        check_cuda_error(cusparseDestroySpMat(persistent_data->matA));
    }
    if(persistent_data->vec_beta_im) {
        check_cuda_error(cusparseDestroyDnVec(persistent_data->vec_beta_im));
    }
    if(persistent_data->vec_vm) {
        check_cuda_error(cusparseDestroyDnVec(persistent_data->vec_vm));
    }
#else
    check_cuda_error((cudaError_t)cusparseDestroyMatDescr(persistent_data->descr));
#endif
    check_cuda_error(cudaFree(persistent_data->d_col));
    check_cuda_error(cudaFree(persistent_data->d_row));
    check_cuda_error(cudaFree(persistent_data->d_val));
    check_cuda_error(cudaFree(persistent_data->beta_im));
    check_cuda_error(cudaFree(persistent_data->tmp_data));
    check_cuda_error(cudaFree(persistent_data->d_distances));
    check_cuda_error(cudaFree(persistent_data->d_volumes));

    fclose(persistent_data->output_file);

    free(persistent_data);
}
#endif

INIT_CALC_ECG(init_pseudo_bidomain) {
    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config, "use_gpu");
    if(gpu) {
#ifdef COMPILE_CUDA
        init_pseudo_bidomain_gpu(config, the_ode_solver, the_grid);
#else
        log_warn("Cuda runtime not found in this system. Falling back to CPU version!!\n");
        init_pseudo_bidomain_cpu(config, NULL, the_grid);
#endif
    } else {
        init_pseudo_bidomain_cpu(config, NULL, the_grid);
    }
}

CALC_ECG(pseudo_bidomain) {
    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config, "use_gpu");
    if(gpu) {
#ifdef COMPILE_CUDA
        pseudo_bidomain_gpu(time_info, config, the_grid);
#else
        pseudo_bidomain_cpu(time_info, config, the_grid);
#endif
    } else {
        pseudo_bidomain_cpu(time_info, config, the_grid);
    }
}

END_CALC_ECG(end_pseudo_bidomain) {
    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config, "use_gpu");
    if(gpu) {
#ifdef COMPILE_CUDA
        end_pseudo_bidomain_gpu(config);
#else
        end_pseudo_bidomain_cpu(config);
#endif
    } else {
        end_pseudo_bidomain_cpu(config);
    }
}
