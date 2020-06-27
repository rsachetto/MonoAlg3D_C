#include "utils.h"

void usage (const char pname[])
{
  fprintf(stderr,"%s\n",PRINT_LINE);
  fprintf(stderr,"Usage:> %s <mode_name>\n",pname);
  fprintf(stderr,"%s\n",PRINT_LINE);
  fprintf(stderr,"<mode_name> = apd | cv | delay\n");
  fprintf(stderr,"\tapd = APD calculator mode (Action Potential Duration)\n");
  fprintf(stderr,"\tcv = CV calculator mode (Conduction Velocity)\n");
  fprintf(stderr,"\tdelay = Delay calculator mode (Time Difference)\n");
  fprintf(stderr,"%s\n",PRINT_LINE);
/*
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
*/
}

void usage_apd (const char pname[])
{
  fprintf(stderr,"%s\n",PRINT_LINE);
  fprintf(stderr,"Usage:> %s apd <grid_filename> <apd_filename>\n",pname);
  fprintf(stderr,"%s\n",PRINT_LINE);
  fprintf(stderr,"<grid_filename> = Filename with the grid configuration\n");
  fprintf(stderr,"<apd_filename> = Filename with the APD of each cell from the grid\n");
  fprintf(stderr,"%s\n",PRINT_LINE);
}

void usage_cv (const char pname[])
{
  fprintf(stderr,"%s\n",PRINT_LINE);
  fprintf(stderr,"Usage:> %s cv <at_filename> <src_index> <dest_index>\n",pname);
  fprintf(stderr,"%s\n",PRINT_LINE);
  fprintf(stderr,"<at_filename> = Filename with the activation time from the grid\n");
  fprintf(stderr,"<src_index> = Cell index of the first reference cell\n");
  fprintf(stderr,"<dest_index> = Cell index of the second reference cell\n");
  fprintf(stderr,"%s\n",PRINT_LINE);
}

std::ostream& operator<<(std::ostream& os, const Point_3D& obj)
{
    os << "(" << obj.x << "," << obj.y << "," << obj.z << ")";
    return os;
}

uint32_t check_input_mode (char *argv[])
{
  char *modename = argv[1];

  if (strcmp(modename,"apd") == 0)
    return 1;
  else if (strcmp(modename,"cv") == 0)
    return 2;
  else if (strcmp(modename,"delay") == 0)
    return 3;
  else
  {
    usage(argv[0]);
    exit(EXIT_FAILURE);
  }
}
