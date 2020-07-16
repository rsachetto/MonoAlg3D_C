// Author: Lucas Berg

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cstdio>

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
#include <vtkLine.h>

#define PRINT_LINE "============================================================================="

using namespace std;

struct cell
{
    uint32_t id;
    double x;
    double y;
    double z;
    double apd;

    bool operator< (const cell &c) const
    {
        return x < c.x;
    }
};


int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        cerr << PRINT_LINE << endl;
        cerr << "Usage:> " << argv[0] << " <apd_filename> <side_length> <dx>" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "<apd_filename> = Filename with the APD from each cell" << endl;
        cerr << "<side_length> = Side length of the grid in {um}" << endl;
        cerr << "<dx> = Mesh discretization in {um}" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " inputs/apd-map.vtu 10000 100" << endl;
        cerr << PRINT_LINE << endl;

        exit(EXIT_FAILURE);
    }

    string filename = argv[1];
    double side_length = atof(argv[2]);
    double dx = atof(argv[3]);

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkUnstructuredGrid *unstructured_grid_2 = reader->GetOutput();

    // Read cells
    uint32_t num_cells = unstructured_grid_2->GetNumberOfCells();
    uint32_t num_points = unstructured_grid_2->GetNumberOfPoints();
    vector<struct cell> the_cells;

    //cout << "Number of points = " << num_points << endl;
    //cout << "Number of cells = " << num_cells << endl;

    for (int i = 0; i < num_cells; i++)
    {
        //cout << "Cell " << i << endl;

        vtkCell *cell = unstructured_grid_2->GetCell(i);

        vtkHexahedron *hexahedron = dynamic_cast<vtkHexahedron*>(cell);

        // Capture the points from the current cell
        double pos_0[3];
        double pos_1[3];
        double pos_3[3];
        double pos_4[3];

        hexahedron->GetPoints()->GetPoint(0,pos_0);
        hexahedron->GetPoints()->GetPoint(1,pos_1);
        hexahedron->GetPoints()->GetPoint(3,pos_3);
        hexahedron->GetPoints()->GetPoint(4,pos_4);

        // Calculates its center 
        double x_bar = (pos_1[0] + pos_0[0]) / 2.0;
        double y_bar = (pos_3[1] + pos_0[1]) / 2.0;
        double z_bar = (pos_4[2] + pos_0[2]) / 2.0;

        struct cell c;
        c.id = i;
        c.x = x_bar;
        c.y = y_bar;
        c.z = z_bar;

        the_cells.push_back(c);
        //cout << "\t\tCentroid " << i << " = (" << x_bar << "," << y_bar << "," << z_bar << ")" << endl;
    }

    // Read cells scalar values
    string array_name = "Scalars_";

    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructured_grid_2->GetCellData()->GetArray(array_name.c_str()));

    if(array)
    {
        // Pass through each cell extracting the scalar value
        for(int i = 0; i < num_cells; i++)
        {
            double value = array->GetValue(i);

            the_cells[i].apd = value;

            //cout << "Cell " << i << " -- Scalar value = " << value << endl;
        }
    }

    // Sort the vector by the x coordinate
    sort(the_cells.begin(),the_cells.end());

    // Write the mean APD per column into a file
    FILE *file = fopen("outputs/mean_apd.txt","w+");
    for (uint32_t i = 0; i < the_cells.size(); i++)
    {
	double x = dx*i + (dx/2.0);
	double apd = the_cells[i].apd;

	fprintf(file,"%g %g\n", x, apd);
    }
    fclose(file);

    return 0;
}
