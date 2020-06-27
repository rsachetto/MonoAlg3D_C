#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <iostream>
#include <string>

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include "../utils/utils.h"

class Config
{
public:
    uint32_t mode_number;

    char *grid_filename;
    char *cells_apd_filename;

    char *activation_time_filename;
    uint32_t src_index;
    uint32_t dest_index;


public:
    Config (uint32_t mode_number, char *grid_filename, char *cells_apd_filename)
    {
      this->mode_number = mode_number;
      this->grid_filename = grid_filename;
      this->cells_apd_filename = cells_apd_filename;
    }

    Config (uint32_t mode_number, char *activation_time_filename, uint32_t src, uint32_t dest)
    {
      this->mode_number = mode_number;
      this->activation_time_filename = activation_time_filename;
      this->src_index = src;
      this->dest_index = dest;
    }

    void print ()
    {
      if (mode_number == 1)
      {
        std::cout << "[APD_Config] Grid filename = " << this->grid_filename << std::endl;
        std::cout << "[APD_Config] Cells APD filename = " << this->cells_apd_filename << std::endl;
      }
      else if (mode_number == 2)
      {
        std::cout << "[CV_Config] Activation time filename = " << this->activation_time_filename << std::endl;
        std::cout << "[CV_Config] Source index = " << this->src_index << std::endl;
        std::cout << "[CV_Config] Destination index = " << this->dest_index << std::endl;
      }
    }
};

uint32_t check_input_mode (char *argv[]);
Config* read_user_input (const uint32_t mode, int argc, char *argv[]);

#endif //MONOALG3D_UTILS_H_H
