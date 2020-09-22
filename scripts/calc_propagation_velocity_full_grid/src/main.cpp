// Author: Lucas Berg
// ---------------------------------------------------------------------------------------------
// Program that writes the activation time and conduction velocity map of a given simulation
// The user needs to pass an activation map as an input parameter.
// This program only works with plain mesh until now ...
// ---------------------------------------------------------------------------------------------

#include <iostream>
#include <vector>

#include "utils/utils.h"
#include "tissue/tissue.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 != 13)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // Input parameters
    string activation_map_filename = argv[1];
    double dt = atof(argv[2]);
    uint32_t print_rate = atoi(argv[3]);
    double side_length_x = atof(argv[4]);
    double side_length_y = atof(argv[5]);
    double dx = atof(argv[6]);
    double dy = atof(argv[7]);
    double min_x = atof(argv[8]);
    double max_x = atof(argv[9]);
    double min_y = atof(argv[10]);
    double max_y = atof(argv[11]);
    double min_z = atof(argv[12]);
    double max_z = atof(argv[13]);
    
    // Create Tissue structure
    Tissue *the_tissue = new Tissue(dt,side_length_x,side_length_y,dx,dy,print_rate,\
                                    min_x,max_x,min_y,max_y,min_z,max_z);

    // Read the data from all the cells inside the input bounds and store into the Tissue structure
    read_data_from_vtu(the_tissue,activation_map_filename);
    //the_tissue->print();
    
    // Calculate propagation velocity using the Finite Difference Method
    calculate_propagation_velocity(the_tissue);
    //the_tissue->print();

    // Write the activation map and conduction velocity map that are inside the bounded region
    //write_scalar_maps_inside_bounds_to_vtu(the_tissue);
    write_scalar_maps_inside_bounds_to_vtu_v2(the_tissue);

    return 0;
}