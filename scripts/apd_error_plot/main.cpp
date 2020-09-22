// Author: Lucas Berg
// Script to write the values from an input file to the cells of the tissue

#include <iostream>
#include <string>
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

#define PRINT_LINE "============================================================================="

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        cerr << PRINT_LINE << endl;
        cerr << "Usage:> " << argv[0] << " <vtu_file> <input_file> <output_filename>" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "<vtu_file> = VTU file with the data from the simulated tissue" << endl;
        cerr << "<input_file> = Input file with the APD error from the simulation" << endl;
        cerr << "<output_filename> = Output filename with the APD error write into the tissue" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " ../../outputs/plain_100_100_100_tentusscher/V_it_0.vtu inputs/Error_1obj.txt outputs/apd_error.vtu" << endl;
        cerr << argv[0] << " ../../outputs/plain_100_100_100_tentusscher/V_it_0.vtu inputs/Error_2obj.txt outputs/apd_error_2.vtu" << endl;
        cerr << PRINT_LINE << endl;
        
        exit(EXIT_FAILURE);
    }

    string vtu_filename = argv[1];
    string error_filename = argv[2];
    string output_filename = argv[3];

    // Read the data from the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(vtu_filename.c_str());
    reader->Update();

    vtkUnstructuredGrid *unstructuredGrid = reader->GetOutput();

    // Read the data from the error file as <vtkFloatArray>
    double value;
    vtkSmartPointer<vtkFloatArray> values = vtkSmartPointer<vtkFloatArray>::New();
    
    FILE *error_file = fopen(error_filename.c_str(),"r");
    while (fscanf(error_file,"%lf",&value) != EOF)
    {
        values->InsertNextValue(value);
    }
    fclose(error_file);

    // Set the error data on the cellData of the vtkUnstructuredGrid
    unstructuredGrid->GetCellData()->SetScalars(values);
    
    // Write the new grid to a VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(output_filename.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();

    return 0;
}