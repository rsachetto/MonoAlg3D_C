#ifndef CALC_PROPAGATION_VELOCITY_TISSUE_H
#define CALC_PROPAGATION_VELOCITY_TISSUE_H

#include <iostream>
#include <string>
#include <map>

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

// Size of the stencil used to calculate the propagation velocity
static const uint32_t OFFSET = 5;

class Tissue
{
public:
    uint32_t num_cells_in_x;
    uint32_t num_cells_in_y;

    double dt;
    double dx;
    double dy;

    uint32_t print_rate;

    double bounds[6];

    uint32_t num_cells;
    Cell *cells;
public:
    Tissue (const double dt, const double side_length_x, const double side_length_y,\
            const double dx, const double dy, const uint32_t print_rate,\
            const double min_x, const double max_x,\
            const double min_y, const double max_y,\
            const double min_z, const double max_z);
    void print ();
    void read_cells_from_vtu (vtkUnstructuredGrid *unstructured_grid);
};

void read_data_from_vtu (Tissue *the_tissue, const std::string filename);
void read_cells_from_vtu (Tissue *the_tissue, vtkUnstructuredGrid *unstructured_grid);

void get_vertex_positions_from_cell(vtkHexahedron *hexahedron,\
    double p0[], double p1[], double p2[], double p3[], double p4[], double p5[], double p6[], double p7[]);
void calculate_cell_middle_position(double center[], const double p0[], const double p1[], const double p3[], const double p4[]);
void calculate_propagation_velocity (Tissue *the_tissue);

void calculate_instantenous_velocity (Tissue *the_tissue, Cell *the_cell,\
                    const uint32_t i, const uint32_t j);
double center_finite_difference (Tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis);
double forward_finite_difference (Tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis);
double backward_finite_difference (Tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis);

void insert_point_into_map (std::map<Point_3D,uint32_t> &points_map,\
                            vtkSmartPointer<vtkPoints> &points, const double p[],\
                            uint32_t *num_points_inside_bounds);

uint32_t get_point_index_from_map (std::map<Point_3D,uint32_t> points_map, const double p[]);

bool check_position (const uint32_t i, const uint32_t j, const uint32_t num_cells_in_x, const uint32_t num_cells_in_y);

void write_scalar_maps_inside_bounds_to_vtu (Tissue *the_tissue);

void write_scalar_maps_inside_bounds_to_vtu_v2 (Tissue *the_tissue);

void insert_point_into_array (vtkSmartPointer<vtkPoints> &points, const double p[]);

// Sorting functions
int sort_by_position (const void *a, const void *b);
int sort_by_index (const void *a, const void *b);

#endif //MONOALG3D_UTILS_H_H