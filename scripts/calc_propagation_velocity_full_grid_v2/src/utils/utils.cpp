#include "utils.h"

void usage (const char pname[])
{
    printf("%s\n",PRINT_LINE);
    printf("Usage:> %s <input_folder> <dt> <print_rate> <side_length_x> <side_length_y> <dx> <dy>\n",pname);
    printf("%s\n",PRINT_LINE);
    printf("<input_folder> = Folder with the simulation data\n");
    printf("<dt> = Size of the time discretization used on the simulation\n");
    printf("<print_rate> = Rate that the output files were written\n");
    printf("<side_length_x> = Size of the tissue in the x axis\n");
    printf("<side_length_y> = Size of the tissue in the y axis\n");
    printf("<dx> = Size of the space discretization in the x axis\n");
    printf("<dy> = Size of the space discretization in the y axis\n");
    printf("%s\n",PRINT_LINE);
}
