#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>

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
#include "../config/config_parser.h"

#include "hexaedron_cell.h"

class Grid
{
public:
  vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkCellArray> cell_array;
  vtkSmartPointer<vtkFloatArray> cell_data;

public:
  Grid (const uint32_t mode, const char filename[], const char extension_name[], const char data_filename[]);
  void read_grid_from_vtu (const char filename[], const char data_filename[]);
  void calculate_propagation_velocity ();
  void write_grid_as_vtu ();
};

Grid* read_grid (Config *the_config);


void insert_point_into_map (std::map<Point_3D,uint32_t> &points_map,\
                            vtkSmartPointer<vtkPoints> &points, const double p[],\
                            uint32_t *num_points_inside_bounds);
void get_extension_name (const char filename[], char extension_name[]);

#endif //MONOALG3D_UTILS_H_H
