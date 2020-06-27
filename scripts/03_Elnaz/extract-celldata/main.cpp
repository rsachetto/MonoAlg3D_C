// Author: Lucas Berg
// A simple program that reads a VTU file and output the CellData into a TXT file.

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
#include <vtkLine.h>

using namespace std;

int main (int argc, char *argv[])
{
  if (argc-1 != 2)
  {
    fprintf(stderr,"Usage:> %s <input_vtu_file> <output_txt_file>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  char *input_vtu_filename = argv[1];
  char *output_txt_filename = argv[2];

  // READER
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(input_vtu_filename);
  reader->Update();

  vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();

  // Read cells
  uint32_t num_cells = unstructured_grid->GetNumberOfCells();
  uint32_t num_points = unstructured_grid->GetNumberOfPoints();

  //cout << "Number of points = " << num_points << endl;
  //cout << "Number of cells = " << num_cells << endl;

  // Read cells scalar values
  string array_name = "Scalars_";
  vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructured_grid->GetCellData()->GetArray(array_name.c_str()));

  // If we found the CellData proceed ...
  if(array)
  {
      // Pass through each cell extracting the scalar value
      vector<double> output_array;
      for(uint32_t i = 0; i < num_cells; i++)
      {
          double value = array->GetValue(i);
          output_array.push_back(value);

          //cout << "Cell " << i << " -- Scalar value = " << value << endl;
      }

      // WRITER
      FILE *file = fopen(output_txt_filename,"w+");
      for(uint32_t i = 0; i < output_array.size(); i++)
      {
        fprintf(file,"%g\n",output_array[i]);
      }
      fclose(file);
  }

  return 0;
}
