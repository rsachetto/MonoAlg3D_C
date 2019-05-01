// Author: Lucas Berg
// WARNING: This code works only for a plain mesh

// TODO: Refactor this code. Divide into modules ...


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

struct cell
{
    uint32_t id;                // Original index

    double x, y, z;             // Coordinates
    double at;                  // Activation time
    double cv;                  // Conduction velocity

    bool in_boundary_1;         // East
    bool in_boundary_2;         // South
    bool in_boundary_3;         // West
    bool in_boundary_4;         // North
};

struct cell* new_cell ()
{
    struct cell *result = (struct cell*)malloc(sizeof(struct cell));

    return result;
}

struct tissue
{
    uint32_t total_num_cells;
    uint32_t num_cells_in_x;
    uint32_t num_cells_in_y;

    double dx;
    double dy;

    struct cell *cells;
};

struct tissue* new_tissue (const uint32_t total_num_cells,\
                        const double side_length_x, const double side_length_y,\
                        const double dx, const double dy)
{
    struct tissue *result = (struct tissue*)malloc(sizeof(struct tissue));

    result->total_num_cells = total_num_cells;
    result->cells = (struct cell*)malloc(sizeof(struct cell)*total_num_cells);

    result->num_cells_in_x = nearbyint(side_length_x / dx);
    result->num_cells_in_y = nearbyint(side_length_y / dy);

    result->dx = dx;
    result->dy = dy;

    return result;
}

void write_conduction_velocity_map_to_vtu (vtkUnstructuredGrid *unstructuredGrid, struct tissue *the_tissue, const string output_filename)
{
    uint32_t total_num_cells = the_tissue->total_num_cells;
    struct cell *cells = the_tissue->cells;

    // Copy the conduction velocity from the tissue structure to a <vtkFloatArray>
    double value;
    vtkSmartPointer<vtkFloatArray> values = vtkSmartPointer<vtkFloatArray>::New();
    
    for (uint32_t i = 0; i < total_num_cells; i++)
    {
        values->InsertNextValue(cells[i].cv);
    }

    // Set the conduction velocity data as the scalars of the CellData from the vtkUnstructuredGrid
    unstructuredGrid->GetCellData()->SetScalars(values);
    
    // Write the new grid to a VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(output_filename.c_str());
    writer->SetInputData(unstructuredGrid);
    writer->Write();

    printf("[+] Conduction velocity map write into '%s' file!\n",output_filename.c_str());
}

void read_data_from_files (struct tissue *the_tissue,\
                            const string position_filename, const string at_filename)
{
    double at;
    double x, y, z;

    FILE *position_file = fopen(position_filename.c_str(),"r");
    FILE *at_file = fopen(at_filename.c_str(),"r");

    uint32_t total_num_cells = the_tissue->total_num_cells;

    for (uint32_t i = 0; i < total_num_cells; i++)
    {
        fscanf(position_file,"%lf %lf %lf",&x,&y,&z);
        fscanf(at_file,"%lf",&at);

        the_tissue->cells[i].id = i;
        the_tissue->cells[i].x = x;
        the_tissue->cells[i].y = y;
        the_tissue->cells[i].z = z;
        the_tissue->cells[i].at = at;

        the_tissue->cells[i].in_boundary_1 = false;
        the_tissue->cells[i].in_boundary_2 = false;
        the_tissue->cells[i].in_boundary_3 = false;
        the_tissue->cells[i].in_boundary_4 = false;
    }

    fclose(at_file);
    fclose(position_file);
}

void print_tissue (struct tissue *the_tissue)
{
    uint32_t nc = the_tissue->total_num_cells;
    struct cell *cells = the_tissue->cells;

    for (uint32_t i = 0; i < nc; i++)
    {
        printf("Cell %u = (%g,%g,%g) -- AT = %g -- Velocity = %g -- (B1 = %u, B2 = %u, B3 = %u, B4 = %u)\n",cells[i].id,\
                cells[i].x,cells[i].y,cells[i].z,cells[i].at,cells[i].cv,\
                cells[i].in_boundary_1,cells[i].in_boundary_2,cells[i].in_boundary_3,cells[i].in_boundary_4);
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

    if (c1->id < c2->id)
        return true;
    else
        return false;
}

void set_boundaries (struct cell *the_cell, const uint32_t i, const uint32_t j,\
                    const uint32_t nx, const uint32_t ny)
{
    // Inside east boundary
    if (i == 0)
        the_cell->in_boundary_1 = true;
    // Inside south boundary
    if (j == (ny-1))
        the_cell->in_boundary_2 = true;
    // Inside west boundary
    if (i == (nx-1))
        the_cell->in_boundary_3 = true;
    // Inside north boundary
    if (j == 0)
        the_cell->in_boundary_4 = true;
}

double center_finite_difference (struct tissue *the_tissue, const uint32_t i, const uint32_t j, const char axis)
{
    double result;

    uint32_t north;
    uint32_t east;
    uint32_t south;
    uint32_t west;

    struct cell *cells = the_tissue->cells;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    if (axis == 'x')
    {
        east = (i-1) * num_cells_in_y + j;
        west = (i+1) * num_cells_in_y + j;

        result = (cells[east].at - cells[west].at) / (2.0*dx);
    }
    else if (axis == 'y')
    {
        north = i * num_cells_in_y + (j-1);
        south = i * num_cells_in_y + (j+1);

        result = (cells[north].at - cells[south].at) / (2.0*dy);
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
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    if (axis == 'x')
    {
        center = i * num_cells_in_y + j;
        west = (i+1) * num_cells_in_y + j;

        result = (cells[west].at - cells[center].at) / (dx);
    }
    else if (axis == 'y')
    {
        center = i * num_cells_in_y + j;
        south = i * num_cells_in_y + (j+1);

        result = (cells[south].at - cells[center].at) / (dy);
    }
    else
    {
        printf("[-] ERROR! On 'forward_finite_difference', invalid axis!\n");
        exit(EXIT_FAILURE);
    }

    return result;
}

double backward_finite_difference (struct tissue *the_tissue, const uint32_t i, const uint32_t j, const char axis)
{
    double result;

    uint32_t center;
    uint32_t north;
    uint32_t east;
    uint32_t south;
    uint32_t west;

    struct cell *cells = the_tissue->cells;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;
    double dx = the_tissue->dx;
    double dy = the_tissue->dy;

    if (axis == 'x')
    {
        center = i * num_cells_in_y + j;
        east = (i-1) * num_cells_in_y + j;

        result = (cells[center].at - cells[east].at) / (dx);
    }
    else if (axis == 'y')
    {
        center = i * num_cells_in_y + j;
        north = i * num_cells_in_y + (j-1);

        result = (cells[center].at - cells[north].at) / (dy);
    }
    else
    {
        printf("[-] ERROR! On 'forward_finite_difference', invalid axis!\n");
        exit(EXIT_FAILURE);
    }

    return result;
}

void calculate_instantenous_velocity (struct tissue *the_tissue, struct cell *the_cell, const uint32_t i, const uint32_t j)
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

    // Calculate the norm of velocity vector and stored in 'cv' variable of the cell
    the_cell->cv = sqrt(pow(vx,2.0) + pow(vy,2.0));

}

void calculate_propagation_velocity (struct tissue *the_tissue)
{
    struct cell *cells = the_tissue->cells;
    uint32_t num_cells_in_x = the_tissue->num_cells_in_x;
    uint32_t num_cells_in_y = the_tissue->num_cells_in_y;

    for (uint32_t i = 0; i < num_cells_in_y; i++)
    {
        for (uint32_t j = 0; j < num_cells_in_x; j++)
        {
            uint32_t k = i * num_cells_in_y + j;

            set_boundaries(&cells[k],i,j,num_cells_in_x,num_cells_in_y);

            calculate_instantenous_velocity(the_tissue,&cells[k],i,j);
        }
    }

    //print_tissue(the_tissue);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 8)
    {
        cerr << PRINT_LINE << endl;
        cerr << "Usage:> " << argv[0] << " <vtu_file> <input_at_file> <input_position_file> <output_filename>";
        cerr << " <side_length_x> <side_length_y> <dx> <dy>" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "<vtu_file> = VTU file with the data from the simulated tissue" << endl;
        cerr << "<input_at_file> = Input file with the activation time of every cell from the simulation" << endl;
        cerr << "<input_position_file> = Input file with the positions of every cell from the simulation" << endl;
        cerr << "<output_filename> = Output filename with the activation time error write into the tissue" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " ../../outputs/plain_100_100_100_fhn/V_it_0.vtu inputs/cells-at.txt inputs/cells-positions.txt outputs/activation_time_map.vtu";
        cerr << " 10000 10000 100 100" << endl;
        cerr << PRINT_LINE << endl;
        
        exit(EXIT_FAILURE);
    }

    string vtu_filename = argv[1];
    string at_filename = argv[2];
    string position_filename = argv[3];
    string output_filename = argv[4];
    double side_length_x = atof(argv[5]);
    double side_length_y = atof(argv[6]);
    double dx = atof(argv[7]);
    double dy = atof(argv[8]);

    vtkUnstructuredGrid *unstructuredGrid;

    // Read the data from the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(vtu_filename.c_str());
    reader->Update();

    // Get a reference to the data as a 'vtkUnstructuredGrid' structure
    unstructuredGrid = reader->GetOutput();
    uint32_t total_num_cells = unstructuredGrid->GetNumberOfCells();

    struct tissue *the_tissue = new_tissue(total_num_cells,side_length_x,side_length_y,dx,dy);
    read_data_from_files(the_tissue,position_filename,at_filename);
    //print_tissue(the_tissue);
    
    qsort(the_tissue->cells,the_tissue->total_num_cells,sizeof(struct cell),sort_by_position);
    //print_tissue(the_tissue);

    calculate_propagation_velocity(the_tissue);
    //print_tissue(the_tissue);

    qsort(the_tissue->cells,the_tissue->total_num_cells,sizeof(struct cell),sort_by_index);

    // Now, the array of cells is in the same order as the cells from the VTU file
    write_conduction_velocity_map_to_vtu (unstructuredGrid,the_tissue,output_filename);

    return 0;
}