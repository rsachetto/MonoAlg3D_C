#ifndef CALC_PROPAGATION_VELOCITY_CONTROL_VOLUME_H
#define CALC_PROPAGATION_VELOCITY_CONTROL_VOLUME_H

#include <iostream>
#include <string>

#include <cstdio>
#include <cstdlib>

#include <vtkXMLUnstructuredGridReader.h>
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

// Control volumes are hexahedrons
//                     p6______p7          
//                     /|     /|           
//                  p2/_|__p3/ |          
//                   |  |___|__|          
//                   | /p5  | /p4           
//                   |/_____|/          
//                  p1      p0
struct control_volume
{
    uint32_t id;

    double p0[3];       // Coordinates of the first point
    double p1[3];       // Coordinates of the second point
    double p2[3];       // Coordinates of the third point
    double p3[3];       // Coordinates of the forth point
    double p4[3];       // Coordinates of the fifth point
    double p5[3];       // Coordinates of the sixty point
    double p6[3];       // Coordinates of the seventy point
    double p7[3];       // Coordinates of the eighty point
};

void read_control_volumes_from_vtu (struct control_volume *volumes,\
                        vtkUnstructuredGrid *unstructured_grid);
void calc_control_volume_middle_position(struct control_volume the_volume, double center[]);



#endif //MONOALG3D_UTILS_H_H