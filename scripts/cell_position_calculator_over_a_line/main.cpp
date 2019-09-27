// Author: Lucas Berg
// Script to write to a file the positions from the center of each control volume

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>

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

#define PRINT_LINE "============================================================================="

using namespace std;

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

struct region
{
    double min_x;       // Minimum value of x
    double max_x;       // Maximum value of x
    double min_y;       // Minimum value of y
    double max_y;       // Maximum value of y
    double min_z;       // Minimum value of z
    double max_z;       // Maximum value of z
};

struct region* new_region (const double min_x, const double max_x,\
                        const double min_y, const double max_y,\
                        const double min_z, const double max_z)
{
    struct region *result = (struct region*)malloc(sizeof(struct region));

    result->min_x = min_x;
    result->max_x = max_x;
    result->min_y = min_y;
    result->max_y = max_y;
    result->min_z = min_z;
    result->max_z = max_z;

    return result;
}

struct tissue
{
    uint32_t num_volumes;

    struct control_volume *volumes;
};

struct tissue* new_tissue ()
{
    struct tissue *t = (struct tissue *)malloc(sizeof(struct tissue));

    t->num_volumes = 0;
    t->volumes = NULL;

    return t;
}

void read_control_volumes_from_vtu (struct tissue *the_tissue, vtkUnstructuredGrid *unstructuredGrid)
{
    uint32_t num_cells = unstructuredGrid->GetNumberOfCells();

    the_tissue->num_volumes = num_cells;
    the_tissue->volumes = (struct control_volume*)malloc(sizeof(struct control_volume)*num_cells);

    for (uint32_t i = 0; i < num_cells; i++)
    {
        vtkCell *cell = unstructuredGrid->GetCell(i);

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

        memcpy(the_tissue->volumes[i].p0,p0,sizeof(double)*3);
        memcpy(the_tissue->volumes[i].p1,p1,sizeof(double)*3);
        memcpy(the_tissue->volumes[i].p2,p2,sizeof(double)*3);
        memcpy(the_tissue->volumes[i].p2,p2,sizeof(double)*3);
        memcpy(the_tissue->volumes[i].p3,p3,sizeof(double)*3);
        memcpy(the_tissue->volumes[i].p4,p4,sizeof(double)*3);
        memcpy(the_tissue->volumes[i].p5,p5,sizeof(double)*3);
        memcpy(the_tissue->volumes[i].p6,p6,sizeof(double)*3);
        memcpy(the_tissue->volumes[i].p7,p7,sizeof(double)*3);
    }    
}

void print_control_volumes (struct tissue *the_tissue)
{
    uint32_t nv = the_tissue->num_volumes;
    struct control_volume *volumes = the_tissue->volumes;

    for (uint32_t i = 0; i < nv; i++)
    {
        printf("%s\n",PRINT_LINE);
        printf("Cell %u\n",i);

        printf("p0: (%g,%g,%g)\n",volumes[i].p0[0],volumes[i].p0[1],volumes[i].p0[2]);
        printf("p1: (%g,%g,%g)\n",volumes[i].p1[0],volumes[i].p1[1],volumes[i].p1[2]);
        printf("p2: (%g,%g,%g)\n",volumes[i].p2[0],volumes[i].p2[1],volumes[i].p2[2]);
        printf("p3: (%g,%g,%g)\n",volumes[i].p3[0],volumes[i].p3[1],volumes[i].p3[2]);
        printf("p4: (%g,%g,%g)\n",volumes[i].p4[0],volumes[i].p4[1],volumes[i].p4[2]);
        printf("p5: (%g,%g,%g)\n",volumes[i].p5[0],volumes[i].p5[1],volumes[i].p5[2]);
        printf("p6: (%g,%g,%g)\n",volumes[i].p6[0],volumes[i].p6[1],volumes[i].p6[2]);
        printf("p7: (%g,%g,%g)\n",volumes[i].p7[0],volumes[i].p7[1],volumes[i].p7[2]);

        printf("%s\n",PRINT_LINE);
    }
}

void calc_control_volume_center (const struct control_volume volume, double center[])
{
    center[0] = (volume.p0[0] + volume.p1[0]) / 2.0;
    center[1] = (volume.p0[1] + volume.p3[1]) / 2.0;
    center[2] = (volume.p0[2] + volume.p4[2]) / 2.0;
}

void write_control_volumes_middle_positions_inside_region (struct tissue *the_tissue, struct region *the_region)
{
    uint32_t nv = the_tissue->num_volumes;
    struct control_volume *volumes = the_tissue->volumes;

    double min_x = the_region->min_x;
    double max_x = the_region->max_x;
    double min_y = the_region->min_y;
    double max_y = the_region->max_y;
    double min_z = the_region->min_z;
    double max_z = the_region->max_z;

    FILE *file = fopen("outputs/cells_positions_inside_region.txt","w+");

    for (uint32_t i = 0; i < nv; i++)
    {
        double center[3];

        calc_control_volume_center(volumes[i],center);

        if (center[0] >= min_x && center[0] <= max_x &&\
            center[1] >= min_y && center[1] <= max_y &&\
            center[2] >= min_z && center[2] <= max_z)
        {
            fprintf(file,"%u %g %g %g\n",i,center[0],center[1],center[2]);
        }
    }

    fclose(file);
}

void create_sphere (vtkSmartPointer<vtkAppendPolyData> appendFilter, const double center[])
{
    // Create the sphere at the specified center
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetCenter(center[0],center[1],center[2]);
    sphereSource->SetRadius(0.05);
    sphereSource->Update();

    // Append the sphere to the filter
    appendFilter->AddInputConnection(sphereSource->GetOutputPort());
    appendFilter->Update();
}

void write_control_volumes_middle_positions_inside_region_to_vtp (struct tissue *the_tissue,\
                                                            struct region *the_region)
{
    uint32_t nv = the_tissue->num_volumes;
    struct control_volume *volumes = the_tissue->volumes;

    double min_x = the_region->min_x;
    double max_x = the_region->max_x;
    double min_y = the_region->min_y;
    double max_y = the_region->max_y;
    double min_z = the_region->min_z;
    double max_z = the_region->max_z;

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    
    // For each control volume create a sphere
    for (uint32_t i = 0; i < nv; i++)
    {
        printf("[main] Working on control volume %u\n",i);

        double center[3];

        calc_control_volume_center(volumes[i],center);

        if (center[0] >= min_x && center[0] <= max_x &&\
            center[1] >= min_y && center[1] <= max_y &&\
            center[2] >= min_z && center[2] <= max_z)
        {
            create_sphere(appendFilter,center);
        }

    }
        
    // Write the data that was appended to the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("outputs/cells_positions_inside_region.vtp");
    writer->SetInputConnection(appendFilter->GetOutputPort());
    writer->Write();
}

void write_control_volumes_cells_indexes_inside_region (struct tissue *the_tissue, struct region *the_region)
{
    uint32_t nv = the_tissue->num_volumes;
    struct control_volume *volumes = the_tissue->volumes;

    double min_x = the_region->min_x;
    double max_x = the_region->max_x;
    double min_y = the_region->min_y;
    double max_y = the_region->max_y;
    double min_z = the_region->min_z;
    double max_z = the_region->max_z;

    FILE *file = fopen("outputs/cells_indexes_inside_region.txt","w+");

    for (uint32_t i = 0; i < nv; i++)
    {
        double center[3];

        calc_control_volume_center(volumes[i],center);

        if (center[0] >= min_x && center[0] <= max_x &&\
            center[1] >= min_y && center[1] <= max_y &&\
            center[2] >= min_z && center[2] <= max_z)
        {
            fprintf(file,"%u\n",i);
        }
    }

    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 7)
    {
        cerr << PRINT_LINE << endl;
        cerr << "Usage:> " << argv[0] << " <vtu_file> <min_x> <max_x> <min_y> <max_y> <min_z> <max_z>" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "<vtu_file> = VTU file with the solution of a timestep from the simulation" << endl;
        cerr << "<min_x> = Minimum value for the x region" << endl;
        cerr << "<max_x> = Maximum value for the x region" << endl;
        cerr << "<min_y> = Minimum value for the y region" << endl;
        cerr << "<max_y> = Maximum value for the y region" << endl;
        cerr << "<min_z> = Minimum value for the z region" << endl;
        cerr << "<max_z> = Maximum value for the z region" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " ../../outputs/plain_mixed_models_tt/V_it_0.vtu" << endl;
        cerr << PRINT_LINE << endl;
        
        exit(EXIT_FAILURE);
    }

    string filename = argv[1];
    double min_x = atof(argv[2]);
    double max_x = atof(argv[3]);
    double min_y = atof(argv[4]);
    double max_y = atof(argv[5]);
    double min_z = atof(argv[6]);
    double max_z = atof(argv[7]);

    // Read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkUnstructuredGrid *unstructuredGrid = reader->GetOutput();

    struct tissue *the_tissue = new_tissue();
    struct region *the_region = new_region(min_x,max_x,min_y,max_y,min_z,max_z);

    read_control_volumes_from_vtu(the_tissue,unstructuredGrid);
    //print_control_volumes(the_tissue);

    write_control_volumes_cells_indexes_inside_region(the_tissue,the_region);
    //write_control_volumes_middle_positions_inside_region(the_tissue,the_region);
    //write_control_volumes_middle_positions_inside_region_to_vtp(the_tissue,the_region);

    return 0;
}
