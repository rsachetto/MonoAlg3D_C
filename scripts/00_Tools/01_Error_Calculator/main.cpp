// Author: Lucas Berg
// Script to write the error values from two input files

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

#define PRINT_LINE "====================================================================================================================================="

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 3)
    {
        cerr << PRINT_LINE << endl;
        cerr << "Usage:> " << argv[0] << " <vtu_file_1> <vtu_file_2> <output_filename>" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "<vtu_file_1> = First VTU file with the data from the simulated tissue" << endl;
        cerr << "<vtu_file_2> = Second VTU file with the data from the simulated tissue" << endl;
        cerr << "<output_filename> = Output filename with data error write into the tissue" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " inputs/activation-map-sc0.vtu inputs/activation-map-sc1.vtu outputs/activation_time_error_s0_s1.vtu" << endl;
        cerr << PRINT_LINE << endl;
        
        exit(EXIT_FAILURE);
    }

    string vtu_filename_1 = argv[1];
    string vtu_filename_2 = argv[2];
    string output_filename = argv[3];

    // Read the data from the first VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader_1 = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader_1->SetFileName(vtu_filename_1.c_str());
    reader_1->Update();

    vtkUnstructuredGrid *unstructured_grid_1 = reader_1->GetOutput();

    // Read the data from the second VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader_2 = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader_2->SetFileName(vtu_filename_2.c_str());
    reader_2->Update();

    vtkUnstructuredGrid *unstructured_grid_2 = reader_2->GetOutput();

    // Get the scalars from the 'vtkUnstructuredGrid' as a 'vtkFloatArray'
    string array_name = "Scalars_";
    vtkIdType total_num_cells_1 = unstructured_grid_1->GetNumberOfCells();
    vtkIdType total_num_cells_2 = unstructured_grid_2->GetNumberOfCells();
    if (total_num_cells_1 != total_num_cells_2)
    {
        cerr << "[-] ERROR! The number of cells are different" << endl;
        exit(EXIT_FAILURE);
    }

    vtkSmartPointer<vtkFloatArray> array_1 = vtkFloatArray::SafeDownCast(unstructured_grid_1->GetCellData()->GetArray(array_name.c_str()));
    vtkSmartPointer<vtkFloatArray> array_2 = vtkFloatArray::SafeDownCast(unstructured_grid_2->GetCellData()->GetArray(array_name.c_str()));

    vtkSmartPointer<vtkFloatArray> values = vtkSmartPointer<vtkFloatArray>::New();

    if(array_1 && array_2)
    {
        // --------------------------------------------------------------------------------------------
        // Calculate global RMS error
        double sum_num = 0.0;
        double sum_den = 0.0;
	
        // Pass through each cell on the tissue for the current timestep
        for(uint32_t i = 0; i < total_num_cells_1; i++)
        {
            double value_1, value_2;
            value_1 = array_1->GetValue(i);
            value_2 = array_2->GetValue(i);

            values->InsertNextValue(fabs(value_1-value_2));

            // RMS calculus
            //double value = (value_1 - value_2) / (value_1);		// Relative error
	        double value = fabs(value_1 - value_2);			        // Absolute error
            sum_num += powf(value,2);                               // RMSE
            sum_den += powf(value_1,2);                             // L2_NORM

        }

        // RMS calculus
        double rmse = sqrt(sum_num/(double)total_num_cells_1);
        double l2_norm = sqrt(sum_den);
        double rrmse = sqrt(sum_num/sum_den);

	//printf("Global RMS = %g\n",rms);
        //printf("RMSE = %g\n",rmse);
        //printf("L2_NORM = %g\n",l2_norm);
        //printf("RRMSE = %g\n",rrmse);
	printf("%g\n",rrmse);
    }
    else
    {
        printf("[-] ERROR! Could not found scalar array '%s' on at least one of the files\n",array_name.c_str());
        exit(EXIT_FAILURE);
    }

    // Set the error data on the cellData of the vtkUnstructuredGrid
    unstructured_grid_1->GetCellData()->SetScalars(values);
    
    // Write the new grid to a VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(output_filename.c_str());
    writer->SetInputData(unstructured_grid_1);
    writer->Write();

    return 0;
}
