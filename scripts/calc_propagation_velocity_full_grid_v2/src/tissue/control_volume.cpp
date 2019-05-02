#include "control_volume.h"

void read_control_volumes_from_vtu (struct control_volume *volumes,\
                        vtkUnstructuredGrid *unstructured_grid)
{
    uint32_t num_cells = unstructured_grid->GetNumberOfCells();

    for (uint32_t i = 0; i < num_cells; i++)
    {
        vtkCell *cell = unstructured_grid->GetCell(i);

        vtkHexahedron *hexahedron = dynamic_cast<vtkHexahedron*>(cell);

        double p0[3];
        double p1[3];
        double p2[3];
        double p3[3];
        double p4[3];
        double p5[3];
        double p6[3];
        double p7[3];

        hexahedron->GetPoints()->GetPoint(0, p0);
        hexahedron->GetPoints()->GetPoint(1, p1);
        hexahedron->GetPoints()->GetPoint(2, p2);
        hexahedron->GetPoints()->GetPoint(3, p3);
        hexahedron->GetPoints()->GetPoint(4, p4);
        hexahedron->GetPoints()->GetPoint(5, p5);
        hexahedron->GetPoints()->GetPoint(6, p6);
        hexahedron->GetPoints()->GetPoint(7, p7);

        memcpy(volumes[i].p0,p0,sizeof(double)*3);
        memcpy(volumes[i].p1,p1,sizeof(double)*3);
        memcpy(volumes[i].p2,p2,sizeof(double)*3);
        memcpy(volumes[i].p2,p2,sizeof(double)*3);
        memcpy(volumes[i].p3,p3,sizeof(double)*3);
        memcpy(volumes[i].p4,p4,sizeof(double)*3);
        memcpy(volumes[i].p5,p5,sizeof(double)*3);
        memcpy(volumes[i].p6,p6,sizeof(double)*3);
        memcpy(volumes[i].p7,p7,sizeof(double)*3);
    }
}

void calc_control_volume_middle_position(struct control_volume the_volume, double center[])
{
    center[0] = (the_volume.p0[0] + the_volume.p1[0]) / 2.0;
    center[1] = (the_volume.p0[1] + the_volume.p3[1]) / 2.0;
    center[2] = (the_volume.p0[2] + the_volume.p4[2]) / 2.0;
}