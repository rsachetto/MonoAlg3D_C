#include "tissue.h"

struct tissue* new_tissue (const double dt, const double side_length_x, const double side_length_y,\
                        const double dx, const double dy, const uint32_t print_rate)
{
    struct tissue *result = (struct tissue*)malloc(sizeof(struct tissue));

    result->num_cells_in_x = nearbyint(side_length_x / dx);
    result->num_cells_in_y = nearbyint(side_length_y / dy);

    result->dt = dt;
    result->dx = dx;
    result->dy = dy;

    result->print_rate = print_rate;

    return result;
}

void set_control_volumes_middle_positions(struct tissue *the_tissue, struct control_volume *volumes)
{
    uint32_t total_num_cells = the_tissue->total_num_cells;
    struct cell *the_cells = the_tissue->cells;

    for (uint32_t i = 0; i < total_num_cells; i++)
    {
        double center[3];

        calc_control_volume_middle_position(volumes[i],center);

        set_new_cell(&the_cells[i],i,center[0],center[1],center[2]);

    }
}

void set_cells_position_with_vtu_file (struct tissue *the_tissue,\
                        const std::string folder_path, const std::string filename)
{
    std::string name = folder_path + "/" + filename;

    // Read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(name.c_str());
    reader->Update();

    vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();

    // Update the total number of cells on the tissue
    the_tissue->total_num_cells = unstructured_grid->GetNumberOfCells();

    // Allocate memory for the cells and control volumes of the tissue
    struct control_volume *volumes = (struct control_volume*)malloc(sizeof(struct control_volume)*the_tissue->total_num_cells);
    the_tissue->cells = (struct cell*)malloc(sizeof(struct cell)*the_tissue->total_num_cells);

    // Read the vertex positions from each control volume on the tissue 
    read_control_volumes_from_vtu(volumes,unstructured_grid);

    // Update the (x,y,z) positions of each cell in the tissue 
    // with middle position from each control volume
    set_control_volumes_middle_positions(the_tissue,volumes);

    free(volumes);

}

void set_activation_times (struct tissue *the_tissue, const std::string folder_path,\
                    std::vector<struct vtu_file> vtu_files)
{   

    uint32_t total_num_cells = the_tissue->total_num_cells;
    uint32_t total_num_steps = vtu_files.size();
    double dt = the_tissue->dt;
    uint32_t print_rate = the_tissue->print_rate;

    struct cell_data *data = new_cell_data(total_num_cells,total_num_steps);

    read_transmembrane_potential_from_vtu(data,folder_path,vtu_files);

    // For each cell calculate the timestep where the maximum derivative occurs
    for (uint32_t i = 0; i < total_num_cells; i++)
    {
        uint32_t max_timestep = -1;
        double max_dvdt = __DBL_MIN__;

        for (uint32_t k = 0; k < total_num_steps-1; k++)
        {
            double dvdt = (data[k+1].vms[i] - data[k].vms[i]) / dt;

            if (dvdt > max_dvdt)
            {
                max_dvdt = dvdt;
                max_timestep = k;
            }
        }

        //printf("Cell %u -- Max timestep = %u\n",i,max_timestep);
        
        // Get the correct time where the maximum derivative occurs for the current cell
        the_tissue->cells[i].at = max_timestep * dt * print_rate;
         
    }

    free_cell_data(data);

    printf("\n[!] Calculated activation time of every cell in the tissue !\n");

}

void load_activation_times (struct tissue *the_tissue, const std::string activation_map_filename)
{
    printf("\n[!] Loading activation map of the tissue from '%s' ...\n",activation_map_filename.c_str());

    std::string array_name = "Scalars_";

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    // Load with a Reader
    reader->SetFileName(activation_map_filename.c_str());
    reader->Update();

    // Get the scalars from the 'vtkUnstructuredGrid' as a 'vtkFloatArray'
    vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();
    vtkIdType total_num_cells = unstructured_grid->GetNumberOfCells();
    vtkSmartPointer<vtkFloatArray> array = 
            vtkFloatArray::SafeDownCast(unstructured_grid->GetCellData()->GetArray(array_name.c_str()));

    if(array)
    {
        // Pass through the transmembrane of each cell on the tissue for the current timestep
        for(int i = 0; i < total_num_cells; i++)
        {
            double value;
            value = array->GetValue(i);

            the_tissue->cells[i].at = value;
            //printf("%u = %g\n",i,value);
        }
    }
    else
    {
        printf("[-] ERROR! Could not found scalar array '%s' on file '%s'\n",array_name.c_str(),activation_map_filename.c_str());
        exit(EXIT_FAILURE);
    }
}

void set_conduction_velocity (struct tissue *the_tissue)
{
    struct cell *cells = the_tissue->cells;
    uint32_t total_num_cells = the_tissue->total_num_cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;

    // Sort the cells by (x, y, z)
    qsort(cells,total_num_cells,sizeof(struct cell),sort_by_position);

    for (uint32_t i = 0; i < num_cells_in_y; i++)
    {
        for (uint32_t j = 0; j < num_cells_in_x; j++)
        {
            uint32_t k = i * num_cells_in_y + j;

            set_boundaries(&cells[k],i,j,num_cells_in_x,num_cells_in_y);

            calculate_instantenous_velocity(the_tissue,&cells[k],i,j);
        }
    }

    // Sort the cells back its original state
    qsort(cells,total_num_cells,sizeof(struct cell),sort_by_index);

    printf("\n[!] Calculated conduction velocity of every cell in the tissue !\n");

}

void read_transmembrane_potential_from_vtu(struct cell_data *the_data, const std::string folder_path,\
                    std::vector<struct vtu_file> vtu_files)
{
    printf("\n[!] Reading transmembrane potential of every cell in the tissue ...\n");

    std::string array_name = "Scalars_";
    uint32_t total_num_steps = the_data->num_steps;

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    // For each timestep read the transmembrane potential of the corresponding 'vtu' file
    for (uint32_t k = 0; k < total_num_steps; k++)
    {
        // Build the name of path_name of the file
        std::string name = folder_path + "/" + vtu_files[k].name;

        // Load with a Reader
        reader->SetFileName(name.c_str());
        reader->Update();

        // Get the scalars from the 'vtkUnstructuredGrid' as a 'vtkFloatArray'
        vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();
        vtkIdType total_num_cells = unstructured_grid->GetNumberOfCells();
        vtkSmartPointer<vtkFloatArray> array = 
                vtkFloatArray::SafeDownCast(unstructured_grid->GetCellData()->GetArray(array_name.c_str()));
    
        if(array)
        {
            // Pass through the transmembrane of each cell on the tissue for the current timestep
            for(int i = 0; i < total_num_cells; i++)
            {
                double value;
                value = array->GetValue(i);

                the_data[k].vms[i] = value;
                //printf("%u = %g\n",i,value);
            }
        }
        else
        {
            printf("[-] ERROR! Could not found scalar array '%s' on file '%s'\n",array_name.c_str(),name.c_str());
            exit(EXIT_FAILURE);
        }
    }

}

void calculate_instantenous_velocity (struct tissue *the_tissue, struct cell *the_cell,\
                    const uint32_t i, const uint32_t j)
{
    double vx, vy;

    bool in_boundary_1 = the_cell->in_boundary_1;
    bool in_boundary_2 = the_cell->in_boundary_2;
    bool in_boundary_3 = the_cell->in_boundary_3;
    bool in_boundary_4 = the_cell->in_boundary_4;

    // Case 1: Interior cell 
    if (!in_boundary_1 && !in_boundary_2 && !in_boundary_3 && !in_boundary_4)
    {
        vx = center_finite_difference(the_tissue,i,j,'x');
        vy = center_finite_difference(the_tissue,i,j,'y');
    }
    // Case 2: Upper right corner
    else if (in_boundary_1 && in_boundary_4)
    {
        vx = forward_finite_difference(the_tissue,i,j,'x');
        vy = forward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 3: Down right corner
    else if (in_boundary_1 && in_boundary_2)
    {
        vx = forward_finite_difference(the_tissue,i,j,'x');
        vy = backward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 4: Down left corner
    else if (in_boundary_2 && in_boundary_3)
    {
        vx = backward_finite_difference(the_tissue,i,j,'x');
        vy = backward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 5: Upper left corner
    else if (in_boundary_3 && in_boundary_4)
    {
        vx = backward_finite_difference(the_tissue,i,j,'x');
        vy = forward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 6: Right boundary
    else if (in_boundary_1)
    {
        vx = forward_finite_difference(the_tissue,i,j,'x');
        vy = center_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 7: Down boundary
    else if (in_boundary_2)
    {
        vx = center_finite_difference(the_tissue,i,j,'x');
        vy = backward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 8: Left boundary
    else if (in_boundary_3)
    {
        vx = backward_finite_difference(the_tissue,i,j,'x');
        vy = center_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 9: Upper boundary
    else if (in_boundary_4)
    {
        vx = center_finite_difference(the_tissue,i,j,'x');
        vy = forward_finite_difference(the_tissue,i,j,'y'); 
    }
    else
    {
        printf("[-] ERROR! On 'calculate_instantenous_velocity', invalid position (i,j)!\n");
        exit(EXIT_FAILURE);
    }

    // The value is given in {ms/um}, 
    // so we need to invert these values to get the correct velocity unit {um/ms}
    vx = 1.0 / vx;
    vy = 1.0 / vy;

    // Check if the value is to big
    if (fabs(vx) > THRESHOLD) vx = 0.0;
    if (fabs(vy) > THRESHOLD) vy = 0.0;

    printf("Cell %u -- v = (%g,%g)\n",the_cell->id,vx,vy);

    // Calculate the norm of velocity vector and stored in 'cv' variable of the cell
    the_cell->cv = sqrt(pow(vx,2.0) + pow(vy,2.0));
}

double center_finite_difference (struct tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis)
{

    struct cell *cells = the_tissue->cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    // Check if the neighbour cells are out of the grid
    //bool north_ok = true;
    //bool east_ok = true;
    //bool south_ok = true;
    //bool west_ok = true;
    bool north_ok = check_position(i,j-OFFSET,num_cells_in_x,num_cells_in_y);
    bool east_ok = check_position(i-OFFSET,j,num_cells_in_x,num_cells_in_y);
    bool south_ok = check_position(i,j+OFFSET,num_cells_in_x,num_cells_in_y);
    bool west_ok = check_position(i+OFFSET,j,num_cells_in_x,num_cells_in_y);

    double result;

    uint32_t north;
    uint32_t east;
    uint32_t south;
    uint32_t west;

    if (axis == 'x')
    {
        east = (i-OFFSET) * num_cells_in_y + j;
        west = (i+OFFSET) * num_cells_in_y + j;

        result = (cells[east].at - cells[west].at) / (2.0*OFFSET*dx);
    }
    else if (axis == 'y')
    {
        north = i * num_cells_in_y + (j-OFFSET);
        south = i * num_cells_in_y + (j+OFFSET);

        result = (cells[north].at - cells[south].at) / (2.0*OFFSET*dy);
    }
    else
    {
        printf("[-] ERROR! On 'center_finite_difference', invalid axis!\n");
        exit(EXIT_FAILURE);
    }

    return result;
}

double forward_finite_difference (struct tissue *the_tissue, const uint32_t i, const uint32_t j, const char axis)
{
    double result;

    uint32_t center;
    uint32_t north;
    uint32_t east;
    uint32_t south;
    uint32_t west;

    struct cell *cells = the_tissue->cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    // Check if the neighbour cells are out of the grid
    //bool north_ok = true;
    //bool east_ok = true;
    //bool south_ok = true;
    //bool west_ok = true;
    bool north_ok = check_position(i,j-OFFSET,num_cells_in_x,num_cells_in_y);
    bool east_ok = check_position(i-OFFSET,j,num_cells_in_x,num_cells_in_y);
    bool south_ok = check_position(i,j+OFFSET,num_cells_in_x,num_cells_in_y);
    bool west_ok = check_position(i+OFFSET,j,num_cells_in_x,num_cells_in_y);

    if (axis == 'x')
    {
        center = i * num_cells_in_y + j;
        west = (i+OFFSET) * num_cells_in_y + j;

        result = (cells[west].at - cells[center].at) / (dx);
    }
    else if (axis == 'y')
    {
        center = i * num_cells_in_y + j;
        south = i * num_cells_in_y + (j+OFFSET);

        result = (cells[south].at - cells[center].at) / (dy);
    }
    else
    {
        printf("[-] ERROR! On 'forward_finite_difference', invalid axis!\n");
        exit(EXIT_FAILURE);
    }

    return result;
}

double backward_finite_difference (struct tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis)
{
    double result;

    uint32_t center;
    uint32_t north;
    uint32_t east;
    uint32_t south;
    uint32_t west;

    struct cell *cells = the_tissue->cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    // Check if the neighbour cells are out of the grid
    //bool north_ok = true;
    //bool east_ok = true;
    //bool south_ok = true;
    //bool west_ok = true;
    bool north_ok = check_position(i,j-OFFSET,num_cells_in_x,num_cells_in_y);
    bool east_ok = check_position(i-OFFSET,j,num_cells_in_x,num_cells_in_y);
    bool south_ok = check_position(i,j+OFFSET,num_cells_in_x,num_cells_in_y);
    bool west_ok = check_position(i+OFFSET,j,num_cells_in_x,num_cells_in_y);

    if (axis == 'x')
    {
        center = i * num_cells_in_y + j;
        east = (i-OFFSET) * num_cells_in_y + j;

        result = (cells[center].at - cells[east].at) / (dx);
    }
    else if (axis == 'y')
    {
        center = i * num_cells_in_y + j;
        north = i * num_cells_in_y + (j-OFFSET);

        result = (cells[center].at - cells[north].at) / (dy);
    }
    else
    {
        printf("[-] ERROR! On 'forward_finite_difference', invalid axis!\n");
        exit(EXIT_FAILURE);
    }

    return result;
}

void write_scalar_map_to_vtu (struct tissue *the_tissue,\
                const std::string folder_path, const std::string filename, const std::string scalar_name)
{

    std::string name = folder_path + "/" + filename;

    // Read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(name.c_str());
    reader->Update();

    vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();

    uint32_t total_num_cells = the_tissue->total_num_cells;
    struct cell *cells = the_tissue->cells;

    // Set the output name
    std::string output_filename;
    if (scalar_name == "activation_time" || scalar_name == "at" || scalar_name == "a")
        output_filename = "outputs/activation_time_map.vtu";
    else if (scalar_name == "conduction_velocity" || scalar_name == "cv" || scalar_name == "c")
        output_filename = "outputs/conduction_velocity_map.vtu";
    else
    {
        printf("[-] ERROR! Invalid scalar_name '%s'\n",scalar_name.c_str());
        exit(EXIT_FAILURE);
    }

    // Copy the conduction velocity from the tissue structure to a <vtkFloatArray>
    double value;
    vtkSmartPointer<vtkFloatArray> values = vtkSmartPointer<vtkFloatArray>::New();
    
    for (uint32_t i = 0; i < total_num_cells; i++)
    {
        if (scalar_name == "activation_time" || scalar_name == "at" || scalar_name == "a")
        {
            values->InsertNextValue(cells[i].at);
        }
        else if (scalar_name == "conduction_velocity" || scalar_name == "cv" || scalar_name == "c")
        {
            values->InsertNextValue(cells[i].cv);
        }
            
        else
        {
            printf("[-] ERROR! Invalid scalar_name '%s'\n",scalar_name.c_str());
            exit(EXIT_FAILURE);
        }
    }

    // Set the conduction velocity data as the scalars of the CellData from the vtkUnstructuredGrid
    unstructured_grid->GetCellData()->SetScalars(values);
    
    // Write the new grid to a VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(output_filename.c_str());
    writer->SetInputData(unstructured_grid);
    int sucess = writer->Write();
    
    if ((scalar_name == "activation_time" || scalar_name == "at" || scalar_name == "a") && (sucess))
        printf("\n[+] Activation time map write into '%s' file!\n",output_filename.c_str());
    else if ((scalar_name == "conduction_velocity" || scalar_name == "cv" || scalar_name == "c") && sucess)
        printf("\n[+] Conduction velocity map write into '%s' file!\n",output_filename.c_str());
    else
    {
        if (!sucess)
            printf("[-] ERROR! Cannot open output file '%s'! Check if the folder 'outputs' exists!\n",output_filename.c_str());
        else
            printf("[-] ERROR! Invalid scalar_name '%s'\n",scalar_name.c_str());
        exit(EXIT_FAILURE);
    }
    
}

int sort_by_position (const void *a, const void *b)
{
    struct cell *c1 = (struct cell *)a;
    struct cell *c2 = (struct cell *)b;

    if (c1->x < c2->x)
        return true;
    else
    {
        if (c1->x == c2->x)
        {
            if (c1->y < c2->y)
                return true;
            else
            {
                if (c1->y == c2->y)
                {
                    if (c1->z < c2->z)
                        return true;
                    else
                        return false;
                }
                else
                {
                    return false;
                }
            }
        }
        else
        {
            return false;
        }
    }
}

int sort_by_index (const void *a, const void *b)
{
    struct cell *c1 = (struct cell *)a;
    struct cell *c2 = (struct cell *)b;

    if (c1->id > c2->id)
        return true;
    else
        return false;
}

bool check_position (const uint32_t i, const uint32_t j, const uint32_t num_cells_in_x, const uint32_t num_cells_in_y)
{
    if (i >= 0 && i < num_cells_in_y && j >= 0 && j < num_cells_in_x)
        return true;
    else
        return false;
}