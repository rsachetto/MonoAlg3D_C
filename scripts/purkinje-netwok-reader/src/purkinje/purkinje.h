#ifndef MONOALG3D_PURKINJE_H
#define MONOALG3D_PURKINJE_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>

// Lines are defined as follows:
//           
//                    ________          
//                  src      dest
struct line
{
    uint32_t id;

    uint32_t src;
    uint32_t dest;
};

struct point
{
    uint32_t id;

    double x, y, z;
};

struct purkinje_network
{
    uint32_t num_points;
    uint32_t num_lines;

    struct point *points;
    struct line *lines;

    double *point_data;
};


struct purkinje_network* new_purkinje_network ();
void read_purkinje_network_from_vtp (struct purkinje_network *the_purkinje_network, vtkPolyData *polydata);
void print_purkinje_network (struct purkinje_network *pk);

struct terminal* calculate_terminals (struct graph *the_graph); 
uint32_t count_number_of_terminals (struct graph *the_graph);

#endif