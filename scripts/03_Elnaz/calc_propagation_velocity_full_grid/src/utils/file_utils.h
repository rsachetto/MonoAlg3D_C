#ifndef CALC_PROPAGATION_VELOCITY_FILE_UTILS_H
#define CALC_PROPAGATION_VELOCITY_FILE_UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <dirent.h>

struct vtu_file
{
    uint32_t timestep;
    std::string name;
};


void read_vtu_files_from_folder(const std::string folder_path,\
                            std::vector<struct vtu_file> &vtu_files);
void filter_files_inside_folder_by_extension (DIR *dir,\
                            const std::string folder_path,const std::string extension_name,\
                            std::vector<struct vtu_file> &vtu_files);
uint32_t get_timestep_from_filename (const std::string filename);
bool sort_by_timestep (struct vtu_file a, struct vtu_file b);

#endif //MONOALG3D_UTILS_H_H