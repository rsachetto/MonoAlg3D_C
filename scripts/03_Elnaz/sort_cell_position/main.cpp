// Author: Lucas Berg

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cstdio>


#define PRINT_LINE "============================================================================="

using namespace std;

struct cell
{
    uint32_t id;
    double x;
    double y;
    double z;

    bool operator< (const cell &c) const
    {
        return x < c.x;
    }
};


int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        cerr << PRINT_LINE << endl;
        cerr << "Usage:> " << argv[0] << " <cell_positions_filename>" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "<cell_positions_filename> = Filename with the cell positions" << endl;
        cerr << PRINT_LINE << endl;
        cerr << "Example:" << endl;
        cerr << argv[0] << " inputs/cell_positions.txt" << endl;
        cerr << PRINT_LINE << endl;
        
        exit(EXIT_FAILURE);
    }

    string filename = argv[1];


    // Read the input file data
    FILE *file = fopen(filename.c_str(),"r");
    uint32_t id;
    double x, y, z;
    vector<struct cell> the_cells;
    while (fscanf(file,"%u %lf %lf %lf",&id,&x,&y,&z) != EOF)
    {
        struct cell c;
        c.id = id;
        c.x = x;
        c.y = y;
        c.z = z;

        the_cells.push_back(c);
    }
    fclose(file);

    // Sort the vector by the x coordinate
    sort(the_cells.begin(),the_cells.end());

    // Write the sorted array on a file
    file = fopen("outputs/sorted_cell_position.txt","w+");
    for (uint32_t i = 0; i < the_cells.size(); i++)
        fprintf(file,"%u %g %g %g\n",the_cells[i].id,the_cells[i].x,the_cells[i].y,the_cells[i].z);
    fclose(file);

    return 0;
}
