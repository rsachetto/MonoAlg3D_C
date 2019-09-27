#include "utils.h"

void usage (const char pname[])
{
    printf("%s\n",PRINT_LINE);
    printf("Usage:> %s <activation_map_filename> <dt> <print_rate> <side_length_x> <side_length_y> <dx> <dy>\n",pname);
    printf("           <min_x> <max_x> <min_y> <max_y> <min_z> <max_z>\n");
    printf("%s\n",PRINT_LINE);
    printf("<activation_map_filename> = Activation map filename\n");
    printf("<dt> = Size of the time discretization used on the simulation\n");
    printf("<print_rate> = Rate that the output files were written\n");
    printf("<side_length_x> = Size of the tissue in the x axis\n");
    printf("<side_length_y> = Size of the tissue in the y axis\n");
    printf("<dx> = Size of the space discretization in the x axis\n");
    printf("<dy> = Size of the space discretization in the y axis\n");
    printf("%s\n",PRINT_LINE);
    printf("Example:\n");
    printf("%s inputs/activation-map.vtu 0.02 50 10000 10000 100 100 4000 8000 3000 7000 0 100\n",pname);
    printf("%s\n",PRINT_LINE);
}


std::ostream& operator<<(std::ostream& os, const Point_3D& obj)
{
    os << "(" << obj.x << "," << obj.y << "," << obj.z << ")";
    return os;
}