// Author: Lucas Berg
// ---------------------------------------------------------------------------------------------
// Program that writes the activation time and conduction velocity map of a given simulation
// The user can pass an activation map as an extra input parameter.
// This program only works with plain mesh until now ...
// ---------------------------------------------------------------------------------------------

#include <iostream>
#include <vector>

#include "utils/utils.h"
#include "tissue/tissue.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 < 7)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    string input_folder = argv[1];
    double dt = atof(argv[2]);
    uint32_t print_rate = atoi(argv[3]);
    double side_length_x = atof(argv[4]);
    double side_length_y = atof(argv[5]);
    double dx = atof(argv[6]);
    double dy = atof(argv[7]);
    
    // Fill the vector with the VTU files from the simulation
    vector<struct vtu_file> vtu_files;
    read_vtu_files_from_folder(input_folder,vtu_files);

    // Configure the tissue structure to store the positions (x,y,z), activation time and
    // conduction velocity for each cell from the tissue
    struct tissue *the_tissue = new_tissue(dt,side_length_x,side_length_y,dx,dy,print_rate);

    // Calculate the cells position (x,y,z)
    set_cells_position_with_vtu_file(the_tissue,input_folder,vtu_files[0].name);

    // If the user passed the file with activation map we will load it
    if (argc-1 == 8)
    {
        string activation_map_filename = argv[8];

        load_activation_times(the_tissue,activation_map_filename);
    }
    // Otherwise, we need to calculate the activation map by hand ...
    else
    {
        set_activation_times(the_tissue,input_folder,vtu_files);
    }

    write_scalar_map_to_vtu(the_tissue,input_folder,vtu_files[0].name,"a");

    set_conduction_velocity(the_tissue);
    write_scalar_map_to_vtu(the_tissue,input_folder,vtu_files[0].name,"c");

    return 0;
}