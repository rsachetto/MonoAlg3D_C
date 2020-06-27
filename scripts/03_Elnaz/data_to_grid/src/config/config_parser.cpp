#include "config_parser.h"

Config* read_user_input (const uint32_t mode, int argc, char *argv[])
{
  Config *result;

  if (mode == 1)
  {
    if (argc-1 != 3)
    {
      usage_apd(argv[0]);
      exit(EXIT_FAILURE);
    }

    result = new Config(mode,argv[2],argv[3]);
    result->print();

  }
  else if (mode == 2)
  {
    if (argc-1 != 4)
    {
      usage_cv(argv[0]);
      exit(EXIT_FAILURE);
    }

    result = new Config(mode,argv[2],atoi(argv[3]),atoi(argv[4]));
    result->print();

  }
  else if (mode == 3)
  {
    // TODO: Need to implement ...
  }
  else
  {
    fprintf(stderr,"[-] ERROR! Invalid mode number '%u'\n",mode);
    exit(EXIT_FAILURE);
  }

  return result;
}
