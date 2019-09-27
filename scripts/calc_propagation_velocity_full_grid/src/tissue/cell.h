#ifndef CALC_PROPAGATION_VELOCITY_CELL_H
#define CALC_PROPAGATION_VELOCITY_CELL_H

#include <iostream>
#include <string>

#include <cstdio>
#include <cstdlib>

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

class Cell
{
public:
    uint32_t id;                // Original index

    double x, y, z;             // Coordinates of the center
    double at;                  // Activation time
    double cv;                  // Conduction velocity

    bool in_boundary_1;         // East
    bool in_boundary_2;         // South
    bool in_boundary_3;         // West
    bool in_boundary_4;         // North

    double p0[3];               // Coordinates of vertex 0
    double p1[3];               // Coordinates of vertex 1
    double p2[3];               // Coordinates of vertex 2
    double p3[3];               // Coordinates of vertex 3
    double p4[3];               // Coordinates of vertex 4
    double p5[3];               // Coordinates of vertex 5
    double p6[3];               // Coordinates of vertex 6
    double p7[3];               // Coordinates of vertex 7   
public:
    Cell ();
    Cell (const uint32_t id, const double x, const double y, const double z);
    Cell (const uint32_t id, const double x, const double y, const double z,\
        const double p0[], const double p1[], const double p2[], const double p3[], const double p4[],\
        const double p5[], const double p6[], const double p7[]);
    void setIndex (const uint32_t id);
    void setCenter (const double center[]);
    void setVertex (const double p0[], const double p1[], const double p2[], const double p3[],\
                    const double p4[], const double p5[], const double p6[], const double p7[]);
    void setActivationTime (const double value);
    void setBoundaries (const uint32_t i, const uint32_t j, const uint32_t nx, const uint32_t ny);
    void print ();
};

void print_cells (const Cell *cells, const uint32_t num_cells);

bool is_inside_bounds (const double center[],\
                    const double min_x, const double max_x,\
                    const double min_y, const double max_y,\
                    const double min_z, const double max_z);


#endif //MONOALG3D_UTILS_H_H