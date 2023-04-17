//
// Created by bergolho on 02/10/20.
//

#ifndef MONOALG3D_SAVE_MESH_HELPERS_H
#define MONOALG3D_SAVE_MESH_HELPERS_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "../3dparty/sds/sds.h"
#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../utils/utils.h"

#include "../libraries_common/common_data_structures.h"
#include "../vtk_utils/vtk_polydata_grid.h"
#include "../vtk_utils/vtk_unstructured_grid.h"

struct common_persistent_data {

    //Ensigth
    uint32_t file_count;
    uint32_t n_digits;

    //Activation times
    struct point_hash_entry *last_time_v;
    struct point_hash_entry *num_activations;
    struct point_hash_entry *cell_was_active;
    struct point_voidp_hash_entry *activation_times;
    struct point_voidp_hash_entry *apds;

    //VTK or VTK
    struct vtk_unstructured_grid *grid;

    bool first_save_call;
    int print_rate;
    int mesh_print_rate;
};

struct save_as_vtp_persistent_data {
    struct vtk_polydata_grid *grid;
    bool first_save_call;
};

struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data {
    struct vtk_unstructured_grid *grid;
    struct vtk_polydata_grid *grid_purkinje;
    bool first_save_call;
};

struct save_coupling_with_activation_times_persistent_data {

    struct vtk_unstructured_grid *tissue_grid;
    struct point_hash_entry *tissue_last_time_v;
    struct point_hash_entry *tissue_num_activations;
    struct point_hash_entry *tissue_cell_was_active;
    struct point_voidp_hash_entry *tissue_activation_times;
    struct point_voidp_hash_entry *tissue_apds;

    struct vtk_polydata_grid *purkinje_grid;
    struct point_hash_entry *purkinje_last_time_v;
    struct point_hash_entry *purkinje_num_activations;
    struct point_hash_entry *purkinje_cell_was_active;
    struct point_voidp_hash_entry *purkinje_activation_times;
    struct point_voidp_hash_entry *purkinje_apds;

    bool first_save_call;
};

struct save_one_cell_state_variables_persistent_data {
    FILE *file;
    char *file_name;
    real_cpu cell_center_x;
    real_cpu cell_center_y;
    real_cpu cell_center_z;
    uint32_t cell_sv_position;
};

struct save_multiple_cell_state_variables_persistent_data {
    uint32_t num_cells;
    FILE **files;
    char *file_name_prefix;
    real_cpu *cell_centers;
    uint32_t *cell_sv_positions;
};

struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data {
    uint32_t num_tissue_cells;
    FILE **tissue_files;
    char *tissue_file_name_prefix;
    real_cpu *tissue_cell_centers;
    uint32_t *tissue_cell_sv_positions;

    uint32_t num_purkinje_cells;
    FILE **purkinje_files;
    char *purkinje_file_name_prefix;
    real_cpu *purkinje_cell_centers;
    uint32_t *purkinje_cell_sv_positions;
};

void add_file_to_pvd(real_cpu current_t, const char *output_dir, const char *base_name, bool first_save_call);
sds create_base_name(char *f_prefix, int iteration_count, char *extension);

// [PURKINJE]
void calculate_purkinje_activation_time_and_apd(struct time_info *time_info, struct config *config, struct grid *the_grid, const real_cpu time_threshold,
                                                const real_cpu purkinje_activation_threshold, const real_cpu purkinje_apd_threshold);
void write_purkinje_activation_time_maps(struct config *config, struct grid *the_grid, char *output_dir, char *file_prefix, bool clip_with_plain,
                                         bool clip_with_bounds, bool binary, bool compress, int compression_level);
void write_purkinje_activation_time_for_each_pulse(struct config *config, struct grid *the_grid, char *output_dir, float plain_coords[], float bounds[],
                                                   bool clip_with_plain, bool clip_with_bounds, bool binary, bool compress, int compression_level);
void write_purkinje_apd_map(struct config *config, struct grid *the_grid, char *output_dir, char *file_prefix, bool clip_with_plain, bool clip_with_bounds,
                            bool binary, bool compress, int compression_level);
void set_purkinje_vtk_values_with_activation_time_from_current_pulse(void **persistent_data, struct grid *the_grid, const int cur_pulse);
void set_purkinje_vtk_values_with_mean_apd(void **persistent_data, struct grid *the_grid);
void print_purkinje_propagation_velocity(struct config *config, struct grid *the_grid);

// [TISSUE]
void calculate_tissue_activation_time_and_apd(struct time_info *time_info, struct config *config, struct grid *the_grid, const real_cpu time_threshold,
                                              const real_cpu tissue_activation_threshold, const real_cpu tissue_apd_threshold);
void write_tissue_activation_time_maps(struct config *config, struct grid *the_grid, char *output_dir, char *file_prefix, bool clip_with_plain,
                                       bool clip_with_bounds, bool binary, bool compress, int compression_level, bool save_f);
void write_tissue_activation_time_for_each_pulse(struct config *config, struct grid *the_grid, char *output_dir, float plain_coords[], float bounds[],
                                                 bool clip_with_plain, bool clip_with_bounds, bool binary, bool compress, int compression_level, bool save_f);
void write_tissue_apd_map(struct config *config, struct grid *the_grid, char *output_dir, char *file_prefix, bool clip_with_plain, bool clip_with_bounds,
                          bool binary, bool compress, int compression_level, bool save_f);
void set_tissue_vtk_values_with_mean_apd(void **persistent_data, struct grid *the_grid);
void set_tissue_vtk_values_with_activation_time_from_current_pulse(void **persistent_data, struct grid *the_grid, const int cur_pulse);

#endif
