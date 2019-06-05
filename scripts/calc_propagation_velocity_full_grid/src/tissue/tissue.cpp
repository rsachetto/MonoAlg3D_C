#include "tissue.h"

void read_data_from_vtu (Tissue *the_tissue, const std::string filename)
{
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();

    // Read the data from every cell in grid and store into 'cells' array 
    the_tissue->read_cells_from_vtu(unstructured_grid);

}

void Tissue::read_cells_from_vtu (vtkUnstructuredGrid *unstructured_grid)
{
    num_cells = unstructured_grid->GetNumberOfCells();
    cells = new Cell[num_cells];

    // Positions from the vertex of each cell
    double p0[3];
    double p1[3];
    double p2[3];
    double p3[3];
    double p4[3];
    double p5[3];
    double p6[3];
    double p7[3];

    // Center position from each cell
    double center[3];

    // Read cells positions
    for (uint32_t i = 0; i < num_cells; i++)
    {
        vtkCell *cell = unstructured_grid->GetCell(i);

        vtkHexahedron *hexahedron = dynamic_cast<vtkHexahedron*>(cell);

        get_vertex_positions_from_cell(hexahedron,p0,p1,p2,p3,p4,p5,p6,p7);

        calculate_cell_middle_position(center,p0,p1,p3,p4);

        // Set the parameters of each Cell in the array
        cells[i].setIndex(i);
        cells[i].setCenter(center);
        cells[i].setVertex(p0,p1,p2,p3,p4,p5,p6,p7);
    }

    // Read cells scalar values
    std::string array_name = "Scalars_";

    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructured_grid->GetCellData()->GetArray(array_name.c_str()));

    if(array)
    {
        // Pass through the activation time of each cell on the tissue
        for(int i = 0; i < num_cells; i++)
        {
            double value = array->GetValue(i);

            cells[i].setActivationTime(value);
        }
    }
    else
    {
        printf("[-] ERROR! Could not found scalar array '%s' on 'vtkUnstructuredGrid'\n",array_name.c_str());
        exit(EXIT_FAILURE);
    }
}

void get_vertex_positions_from_cell(vtkHexahedron *hexahedron,\
    double p0[], double p1[], double p2[], double p3[], double p4[], double p5[], double p6[], double p7[])
{
    hexahedron->GetPoints()->GetPoint(0, p0);
    hexahedron->GetPoints()->GetPoint(1, p1);
    hexahedron->GetPoints()->GetPoint(2, p2);
    hexahedron->GetPoints()->GetPoint(3, p3);
    hexahedron->GetPoints()->GetPoint(4, p4);
    hexahedron->GetPoints()->GetPoint(5, p5);
    hexahedron->GetPoints()->GetPoint(6, p6);
    hexahedron->GetPoints()->GetPoint(7, p7);
}

void calculate_cell_middle_position(double center[], const double p0[], const double p1[], const double p3[], const double p4[])
{
    center[0] = (p0[0] + p1[0]) / 2.0;
    center[1] = (p0[1] + p3[1]) / 2.0;
    center[2] = (p0[2] + p4[2]) / 2.0;
}

void calculate_propagation_velocity (Tissue *the_tissue)
{
    Cell *cells = the_tissue->cells;
    uint32_t num_cells = the_tissue->num_cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;

    // Sort the cells by (x, y, z)
    qsort(cells,num_cells,sizeof(Cell),sort_by_position);
    //print_cells(cells,num_cells);

    for (uint32_t i = 0; i < num_cells_in_y; i++)
    {
        for (uint32_t j = 0; j < num_cells_in_x; j++)
        {
            uint32_t k = i * num_cells_in_y + j;

            cells[k].setBoundaries(i,j,num_cells_in_x,num_cells_in_y);

            calculate_instantenous_velocity(the_tissue,&cells[k],i,j);
        }
    }

    // Sort the cells back its original state
    qsort(cells,num_cells,sizeof(Cell),sort_by_index);

    printf("[!] Calculated conduction velocity of every cell in the tissue !\n");
}

Tissue::Tissue (const double dt, const double side_length_x, const double side_length_y,\
            const double dx, const double dy, const uint32_t print_rate,\
            const double min_x, const double max_x,\
            const double min_y, const double max_y,\
            const double min_z, const double max_z)
{
    this->num_cells_in_x = nearbyint(side_length_x / dx);
    this->num_cells_in_y = nearbyint(side_length_y / dy);

    this->dt = dt;
    this->dx = dx;
    this->dy = dy;

    this->bounds[0] = min_x;
    this->bounds[1] = max_x;
    this->bounds[2] = min_y;
    this->bounds[3] = max_y;
    this->bounds[4] = min_z;
    this->bounds[5] = max_z;

    this->print_rate = print_rate;

}

void Tissue::print ()
{
    printf("======================== TISSUE ============================\n");
    printf("Total number of cells = %u\n",num_cells);
    printf("Number of cells in x = %u\n",num_cells_in_x);
    printf("Number of cells in y = %u\n",num_cells_in_y);
    printf("dt = %g\n",dt);
    printf("dx = %g\n",dx);
    printf("dy = %g\n",dy);
    printf("print_rate = %u\n",print_rate);
    printf("****** BOUNDS ******\n");
    printf("min_x = %g\n",bounds[0]);
    printf("max_x = %g\n",bounds[1]);
    printf("min_y = %g\n",bounds[2]);
    printf("max_y = %g\n",bounds[3]);
    printf("min_z = %g\n",bounds[4]);
    printf("max_z = %g\n",bounds[5]);
    printf("****** CELLS ******\n");
    for (uint32_t i = 0; i < num_cells; i++)
    {
        cells[i].print();
    }

    printf("======================== TISSUE ============================\n");
}

void calculate_instantenous_velocity (Tissue *the_tissue, Cell *the_cell,\
                    const uint32_t i, const uint32_t j)
{
    // Gradient of the activation time
    double grad_at[2];

    bool in_boundary_1 = the_cell->in_boundary_1;
    bool in_boundary_2 = the_cell->in_boundary_2;
    bool in_boundary_3 = the_cell->in_boundary_3;
    bool in_boundary_4 = the_cell->in_boundary_4;

    // Case 1: Interior cell 
    if (!in_boundary_1 && !in_boundary_2 && !in_boundary_3 && !in_boundary_4)
    {
        grad_at[0] = center_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = center_finite_difference(the_tissue,i,j,'y');
    }
    // Case 2: Upper right corner
    else if (in_boundary_1 && in_boundary_4)
    {
        grad_at[0] = forward_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = forward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 3: Down right corner
    else if (in_boundary_1 && in_boundary_2)
    {
        grad_at[0] = forward_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = backward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 4: Down left corner
    else if (in_boundary_2 && in_boundary_3)
    {
        grad_at[0] = backward_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = backward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 5: Upper left corner
    else if (in_boundary_3 && in_boundary_4)
    {
        grad_at[0] = backward_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = forward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 6: Right boundary
    else if (in_boundary_1)
    {
        grad_at[0] = forward_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = center_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 7: Down boundary
    else if (in_boundary_2)
    {
        grad_at[0] = center_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = backward_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 8: Left boundary
    else if (in_boundary_3)
    {
        grad_at[0] = backward_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = center_finite_difference(the_tissue,i,j,'y'); 
    }
    // Case 9: Upper boundary
    else if (in_boundary_4)
    {
        grad_at[0] = center_finite_difference(the_tissue,i,j,'x');
        grad_at[1] = forward_finite_difference(the_tissue,i,j,'y'); 
    }
    else
    {
        printf("[-] ERROR! On 'calculate_instantenous_velocity', invalid position (i,j)!\n");
        exit(EXIT_FAILURE);
    }

    // Calculate the norm of the gradient vector
    double norm_grad_at = sqrt(pow(grad_at[0],2.0) + pow(grad_at[1],2.0));

    // The value is given in {ms/um}, 
    // so we need to invert the value to get the correct velocity unit {um/ms}
    the_cell->cv = pow(norm_grad_at,-1.0);

}

double center_finite_difference (Tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis)
{

    Cell *cells = the_tissue->cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    // Check if the neighbour cells are out of the grid
    bool north_ok = check_position(i,j-OFFSET,num_cells_in_x,num_cells_in_y);
    bool east_ok = check_position(i-OFFSET,j,num_cells_in_x,num_cells_in_y);
    bool south_ok = check_position(i,j+OFFSET,num_cells_in_x,num_cells_in_y);
    bool west_ok = check_position(i+OFFSET,j,num_cells_in_x,num_cells_in_y);

    double result;

    uint32_t center;
    uint32_t north;
    uint32_t east;
    uint32_t south;
    uint32_t west;

    center = i * num_cells_in_y + j;
    east = (i-OFFSET) * num_cells_in_y + j;
    west = (i+OFFSET) * num_cells_in_y + j;
    north = i * num_cells_in_y + (j-OFFSET);
    south = i * num_cells_in_y + (j+OFFSET);

    if (north_ok && east_ok && south_ok && west_ok)
    {
        if (axis == 'x')
        {
            result = (cells[east].at - cells[west].at) / (2.0*OFFSET*dx);
        }
        else if (axis == 'y')
        {
            result = (cells[north].at - cells[south].at) / (2.0*OFFSET*dy);
        }
        else
        {
            printf("[-] ERROR! On 'center_finite_difference', invalid axis!\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        result = -1.0;
    }
    
    return result;
}

double forward_finite_difference (Tissue *the_tissue,\
                            const uint32_t i, const uint32_t j, const char axis)
{
    Cell *cells = the_tissue->cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    // Check if the neighbour cells are out of the grid
    bool north_ok = check_position(i,j-OFFSET,num_cells_in_x,num_cells_in_y);
    bool east_ok = check_position(i-OFFSET,j,num_cells_in_x,num_cells_in_y);
    bool south_ok = check_position(i,j+OFFSET,num_cells_in_x,num_cells_in_y);
    bool west_ok = check_position(i+OFFSET,j,num_cells_in_x,num_cells_in_y);

    double result;

    uint32_t center;
    uint32_t north;
    uint32_t east;
    uint32_t south;
    uint32_t west;

    center = i * num_cells_in_y + j;
    east = (i-OFFSET) * num_cells_in_y + j;
    west = (i+OFFSET) * num_cells_in_y + j;
    north = i * num_cells_in_y + (j-OFFSET);
    south = i * num_cells_in_y + (j+OFFSET);

    if (north_ok && east_ok && south_ok && west_ok)
    {
        if (axis == 'x')
        {
            result = (cells[west].at - cells[center].at) / (dx);
        }
        else if (axis == 'y')
        {
            result = (cells[south].at - cells[center].at) / (dy);
        }
        else
        {
            printf("[-] ERROR! On 'forward_finite_difference', invalid axis!\n");
            exit(EXIT_FAILURE);
        }
    }  
    else
    {
        result = -1.0;
    }
    
    return result;
}

double backward_finite_difference (Tissue *the_tissue,\
                    const uint32_t i, const uint32_t j, const char axis)
{

    Cell *cells = the_tissue->cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    // Check if the neighbour cells are out of the grid
    bool north_ok = check_position(i,j-OFFSET,num_cells_in_x,num_cells_in_y);
    bool east_ok = check_position(i-OFFSET,j,num_cells_in_x,num_cells_in_y);
    bool south_ok = check_position(i,j+OFFSET,num_cells_in_x,num_cells_in_y);
    bool west_ok = check_position(i+OFFSET,j,num_cells_in_x,num_cells_in_y);

    uint32_t center;
    uint32_t north;
    uint32_t east;
    uint32_t south;
    uint32_t west;

    double result;

    if (north_ok && east_ok && south_ok && west_ok)
    {
        if (axis == 'x')
        {
            result = (cells[center].at - cells[east].at) / (dx);
        }
        else if (axis == 'y')
        {
            result = (cells[center].at - cells[north].at) / (dy);
        }
        else
        {
            printf("[-] ERROR! On 'forward_finite_difference', invalid axis!\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        result = -1.0;
    }

    return result;
}

// This version uses a unique map to store the points
void write_scalar_maps_inside_bounds_to_vtu (Tissue *the_tissue)
{
    Cell *cells = the_tissue->cells;
    uint32_t num_cells = the_tissue->num_cells;

    double center[3];

    double min_x = the_tissue->bounds[0];
    double max_x = the_tissue->bounds[1];
    double min_y = the_tissue->bounds[2];
    double max_y = the_tissue->bounds[3];
    double min_z = the_tissue->bounds[4];
    double max_z = the_tissue->bounds[5];

    // Map to store the numbering of the points in a unique way
    std::map<Point_3D,uint32_t> points_map;

    // Array to store all the points inside the bounds
    uint32_t num_points_inside_bounds = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Array to store all the cells inside the bounds
    vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();

    // Arrays to store the activation time and conduction velocity values
    vtkSmartPointer<vtkFloatArray> at_values = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray> cv_values = vtkSmartPointer<vtkFloatArray>::New();

    // For each cell check which ones are inside the bounded region
    for (uint32_t i = 0; i < num_cells; i++)
    {
        center[0] = cells[i].x;
        center[1] = cells[i].y;
        center[2] = cells[i].z;

        if (is_inside_bounds(center,min_x,max_x,min_y,max_y,min_z,max_z))
        {
            // Get a reference for each vertex from the cell
            double *p0 = cells[i].p0;
            double *p1 = cells[i].p1;
            double *p2 = cells[i].p2;
            double *p3 = cells[i].p3;
            double *p4 = cells[i].p4;
            double *p5 = cells[i].p5;
            double *p6 = cells[i].p6;
            double *p7 = cells[i].p7;

            // Try to insert each point from the cell into the map
            insert_point_into_map(points_map,points,p0,&num_points_inside_bounds);
            insert_point_into_map(points_map,points,p1,&num_points_inside_bounds);
            insert_point_into_map(points_map,points,p2,&num_points_inside_bounds);
            insert_point_into_map(points_map,points,p3,&num_points_inside_bounds);
            insert_point_into_map(points_map,points,p4,&num_points_inside_bounds);
            insert_point_into_map(points_map,points,p5,&num_points_inside_bounds);
            insert_point_into_map(points_map,points,p6,&num_points_inside_bounds);
            insert_point_into_map(points_map,points,p7,&num_points_inside_bounds);

            uint32_t id0;
            uint32_t id1;
            uint32_t id2;
            uint32_t id3;
            uint32_t id4;
            uint32_t id5;
            uint32_t id6;
            uint32_t id7;

            // Get the index from the points of each vertex from the cell
            id0 = get_point_index_from_map(points_map,p0);
            id1 = get_point_index_from_map(points_map,p1);
            id2 = get_point_index_from_map(points_map,p2);
            id3 = get_point_index_from_map(points_map,p4);
            id4 = get_point_index_from_map(points_map,p4);
            id5 = get_point_index_from_map(points_map,p5);
            id6 = get_point_index_from_map(points_map,p6);
            id7 = get_point_index_from_map(points_map,p7);

            // Build a Hexahedron cell
            vtkSmartPointer<vtkHexahedron> hexahedron = vtkSmartPointer<vtkHexahedron>::New();
            hexahedron->GetPointIds()->SetId(0, id0);
            hexahedron->GetPointIds()->SetId(1, id1);
            hexahedron->GetPointIds()->SetId(2, id2);
            hexahedron->GetPointIds()->SetId(3, id3);
            hexahedron->GetPointIds()->SetId(4, id4);
            hexahedron->GetPointIds()->SetId(5, id5);
            hexahedron->GetPointIds()->SetId(6, id6);
            hexahedron->GetPointIds()->SetId(7, id7);

            // Add the cell on the cell array
            cell_array->InsertNextCell(hexahedron);

            // Insert the cell values into the arrays
            at_values->InsertNextValue(cells[i].at);
            cv_values->InsertNextValue(cells[i].cv);
        }
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid_new = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructured_grid_new->SetPoints(points);
    unstructured_grid_new->SetCells(VTK_HEXAHEDRON, cell_array);
    unstructured_grid_new->GetCellData()->SetScalars(at_values);

    printf("[+] Writing clipped activation map on '%s'\n","outputs/activation_time_map.vtu");
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("outputs/activation_time_map.vtu");
    writer->SetInputData(unstructured_grid_new);
    writer->Write();

    printf("[+] Writing clipped conduction velocity map on '%s'\n","outputs/conduction_velocity_map.vtu");
    unstructured_grid_new->GetCellData()->SetScalars(cv_values);
    writer->SetFileName("outputs/conduction_velocity_map.vtu");
    writer->Write();
}

// This version adds all the points inside the array without checking if it is already there
void write_scalar_maps_inside_bounds_to_vtu_v2 (Tissue *the_tissue)
{
    Cell *cells = the_tissue->cells;
    uint32_t num_cells = the_tissue->num_cells;

    double center[3];

    double min_x = the_tissue->bounds[0];
    double max_x = the_tissue->bounds[1];
    double min_y = the_tissue->bounds[2];
    double max_y = the_tissue->bounds[3];
    double min_z = the_tissue->bounds[4];
    double max_z = the_tissue->bounds[5];

    // Array to store all the points inside the bounds
    uint32_t num_points_inside_bounds = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Array to store all the cells inside the bounds
    vtkSmartPointer<vtkCellArray> cell_array = vtkSmartPointer<vtkCellArray>::New();

    // Arrays to store the activation time and conduction velocity values
    vtkSmartPointer<vtkFloatArray> at_values = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray> cv_values = vtkSmartPointer<vtkFloatArray>::New();

    // For each cell check which ones are inside the bounded region
    for (uint32_t i = 0; i < num_cells; i++)
    {
        center[0] = cells[i].x;
        center[1] = cells[i].y;
        center[2] = cells[i].z;

        if (is_inside_bounds(center,min_x,max_x,min_y,max_y,min_z,max_z))
        {
            // Get a reference for each vertex from the cell
            double *p0 = cells[i].p0;
            double *p1 = cells[i].p1;
            double *p2 = cells[i].p2;
            double *p3 = cells[i].p3;
            double *p4 = cells[i].p4;
            double *p5 = cells[i].p5;
            double *p6 = cells[i].p6;
            double *p7 = cells[i].p7;

            // Insert each point from the cell into the points array
            insert_point_into_array(points,p0);
            insert_point_into_array(points,p1);
            insert_point_into_array(points,p2);
            insert_point_into_array(points,p3);
            insert_point_into_array(points,p4);
            insert_point_into_array(points,p5);
            insert_point_into_array(points,p6);
            insert_point_into_array(points,p7);

            // Build a Hexahedron cell
            vtkSmartPointer<vtkHexahedron> hexahedron = vtkSmartPointer<vtkHexahedron>::New();
            hexahedron->GetPointIds()->SetId(0, num_points_inside_bounds);
            hexahedron->GetPointIds()->SetId(1, num_points_inside_bounds+1);
            hexahedron->GetPointIds()->SetId(2, num_points_inside_bounds+2);
            hexahedron->GetPointIds()->SetId(3, num_points_inside_bounds+3);
            hexahedron->GetPointIds()->SetId(4, num_points_inside_bounds+4);
            hexahedron->GetPointIds()->SetId(5, num_points_inside_bounds+5);
            hexahedron->GetPointIds()->SetId(6, num_points_inside_bounds+6);
            hexahedron->GetPointIds()->SetId(7, num_points_inside_bounds+7);

            // Add the cell on the cell array
            cell_array->InsertNextCell(hexahedron);

            // Insert the cell values into the arrays
            at_values->InsertNextValue(cells[i].at);
            cv_values->InsertNextValue(cells[i].cv);

            num_points_inside_bounds += 8;
        }
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid_new = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructured_grid_new->SetPoints(points);
    unstructured_grid_new->SetCells(VTK_HEXAHEDRON, cell_array);
    unstructured_grid_new->GetCellData()->SetScalars(at_values);

    printf("[+] Writing clipped activation map on '%s'\n","outputs/activation_time_map.vtu");
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("outputs/activation_time_map.vtu");
    writer->SetInputData(unstructured_grid_new);
    writer->Write();

    printf("[+] Writing clipped conduction velocity map on '%s'\n","outputs/conduction_velocity_map.vtu");
    unstructured_grid_new->GetCellData()->SetScalars(cv_values);
    writer->SetFileName("outputs/conduction_velocity_map.vtu");
    writer->Write();
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

void insert_point_into_array (vtkSmartPointer<vtkPoints> &points, const double p[])
{
    Point_3D point(p[0],p[1],p[2]);

    points->InsertNextPoint(p[0],p[1],p[2]);
}

uint32_t get_point_index_from_map (std::map<Point_3D,uint32_t> points_map, const double p[])
{
    std::map<Point_3D,uint32_t>::iterator it;
    
    Point_3D point(p[0],p[1],p[2]);

    it = points_map.find(point); 
    if (it != points_map.end())
        return it->second;
    else
    {
        printf("[-] ERROR! Point (%g,%g,%g) not found in the points map!\n",p[0],p[1],p[2]);
        exit(EXIT_FAILURE);
    }
}

int sort_by_position (const void *a, const void *b)
{
    Cell *c1 = (Cell*)a;
    Cell *c2 = (Cell*)b;

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
    Cell *c1 = (Cell*)a;
    Cell *c2 = (Cell*)b;

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