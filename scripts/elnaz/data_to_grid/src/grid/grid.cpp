#include "grid.h"

Grid* read_grid (Config *the_config)
{
  char *filename = NULL;
  char *data_filename = NULL;
  uint32_t src_index, dest_index;
  char extension_name[4];
  uint32_t mode = the_config->mode_number;


  // APD: Read the grid configuration
  if (mode == 1)
  {
    filename = the_config->grid_filename;
    data_filename = the_config->cells_apd_filename;
  }
  // CV: Read the grid configuration from the activation time map
  else if (mode == 2)
  {
    filename = the_config->activation_time_filename;

  }
  else if (mode == 3)
  {
    // TODO
  }

  get_extension_name(filename,extension_name);

  Grid *result = new Grid(mode,filename,extension_name,data_filename);

  return result;
}

void get_extension_name (const char filename[], char extension_name[])
{
  uint32_t size = strlen(filename);

  extension_name[0] = filename[size-3];
  extension_name[1] = filename[size-2];
  extension_name[2] = filename[size-1];
  extension_name[4] = '\0';
}


Grid::Grid (const uint32_t mode, const char filename[], const char extension_name[], const char data_filename[])
{
  if (strcmp(extension_name,"vtu") == 0)
  {
    // Data structures to store the grid data
    this->unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    this->points = vtkSmartPointer<vtkPoints>::New();
    this->cell_array = vtkSmartPointer<vtkCellArray>::New();
    this->cell_data = vtkSmartPointer<vtkFloatArray>::New();

    read_grid_from_vtu(filename,data_filename);

    if (mode == 1)
      write_grid_as_vtu();
    else if (mode == 2)
      calculate_propagation_velocity();
  }
  else if (strcmp(extension_name,"vtp") == 0)
  {

  }
}


void Grid::read_grid_from_vtu (const char filename[], const char data_filename[])
{
  // Open the grid file
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename);
  reader->Update();

  vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();

  uint32_t num_points = unstructured_grid->GetNumberOfPoints();
  uint32_t num_cells = unstructured_grid->GetNumberOfCells();

  std::map<Point_3D,uint32_t> points_map;
  uint32_t cur_point = 0;

  // Read all cells
  for (int i = 0; i < num_cells; i++)
  {

    vtkCell *cell = unstructured_grid->GetCell(i);
    vtkHexahedron *hexahedron = dynamic_cast<vtkHexahedron*>(cell);

    // Read points from each cell and store them on a map to avoid duplicates
    for (int j = 0; j < 8; j++)
    {
      double pos[3];
      hexahedron->GetPoints()->GetPoint(j,pos);

      insert_point_into_map(points_map,this->points,pos,&cur_point);
    }

    // Store the current cell on an array
    this->cell_array->InsertNextCell(hexahedron);
  }

  // Read the APD data
  if (data_filename)
  {
    FILE *file = fopen(data_filename,"r");

    double value;
    while (fscanf(file,"%lf",&value) != EOF)
    {
      this->cell_data->InsertNextValue(value);
    }

    fclose(file);
  }
  // Else read the Activation time data
  else
  {
    std::string array_name = "Scalars_";

    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructured_grid->GetCellData()->GetArray(array_name.c_str()));

    if(array)
    {
        // Pass through the activation time of each cell on the tissue
        for(int i = 0; i < num_cells; i++)
        {
            double value = array->GetValue(i);

            this->cell_data->InsertNextValue(value);
        }
    }
    else
    {
        printf("[-] ERROR! Could not found scalar array '%s' on 'vtkUnstructuredGrid'\n",array_name.c_str());
        exit(EXIT_FAILURE);
    }
  }

}

void Grid::calculate_propagation_velocity ()
{

}

void insert_point_into_map (std::map<Point_3D,uint32_t> &points_map,\
                            vtkSmartPointer<vtkPoints> &points, const double p[],\
                            uint32_t *num_points_inside_bounds)
{
    Point_3D point(p[0],p[1],p[2]);

    if (points_map.find(point) == points_map.end())
    {
        points_map.insert(std::pair<Point_3D,uint32_t>(point,*num_points_inside_bounds));

        points->InsertNextPoint(p[0],p[1],p[2]);

        *num_points_inside_bounds = *num_points_inside_bounds + 1;
    }
}

void Grid::write_grid_as_vtu ()
{
  this->unstructured_grid->SetPoints(this->points);
  this->unstructured_grid->SetCells(VTK_HEXAHEDRON,this->cell_array);
  this->unstructured_grid->GetCellData()->SetScalars(this->cell_data);

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName("outputs/tissue_apd_map.vtu");
  writer->SetInputData(this->unstructured_grid);
  writer->Write();
}
