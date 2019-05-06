#ifndef CALC_PROPAGATION_VELOCITY_TISSUE_H
#define CALC_PROPAGATION_VELOCITY_TISSUE_H

#include <iostream>
#include <string>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkHexahedron.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>

#include "../utils/utils.h"

#include "cell.h"
#include "control_volume.h"
#include "cell_data.h"

// Size of the stencil used to calculate the propagation velocity
static const uint32_t OFFSET = 1;

struct tissue
{
    std::string *folder_path;

    uint32_t num_cells_in_x;
    uint32_t num_cells_in_y;

    double dt;
    double dx;
    double dy;

    uint32_t print_rate;

    uint32_t total_num_cells;
    struct cell *cells;
};

struct tissue* new_tissue (const double dt, const double side_length_x, const double side_length_y,\
                        const double dx, const double dy, const uint print_rate);


void set_cells_position_with_vtu_file (struct tissue *the_tissue,\
                        const std::string folder_path, const std::string filename);
void set_control_volumes_middle_positions(struct tissue *the_tissue, struct control_volume *volumes);
void set_activation_times (struct tissue *the_tissue, const std::string folder_path,\
                    std::vector<struct vtu_file> vtu_files);
void set_conduction_velocity (struct tissue *the_tissue);

void load_activation_times (struct tissue *the_tissue, const std::string activation_map_filename);

void read_transmembrane_potential_from_vtu(struct cell_data *the_data, const std::string folder_path,\
                    std::vector<struct vtu_file> vtu_files);


void calculate_instantenous_velocity (struct tissue *the_tissue, struct cell *the_cell,\
                    const uint32_t i, const uint32_t j);
double center_finite_difference (struct tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis);
double forward_finite_difference (struct tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis);
double backward_finite_difference (struct tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis);


void write_scalar_map_to_vtu (struct tissue *the_tissue,\
                const std::string folder_path, const std::string filename, const std::string scalar_name);

int sort_by_position (const void *a, const void *b);
int sort_by_index (const void *a, const void *b);

bool check_position (const uint32_t i, const uint32_t j, const uint32_t num_cells_in_x, const uint32_t num_cells_in_y);

#endif //MONOALG3D_UTILS_H_H