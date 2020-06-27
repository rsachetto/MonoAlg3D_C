#include <iostream>

#include "config/config_parser.h"
#include "utils/utils.h"
#include "grid/grid.h"

int main (int argc, char *argv[])
{
  if (argc-1 < 1)
  {
    usage(argv[0]);
    exit(EXIT_FAILURE);
  }

  uint32_t mode = check_input_mode(argv);

  Config *the_config = read_user_input(mode,argc,argv);

  Grid *the_grid = read_grid(the_config);

}
