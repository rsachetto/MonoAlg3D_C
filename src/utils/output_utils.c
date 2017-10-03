//
// Created by sachetto on 03/10/17.
//

#include <sys/stat.h>
#include "output_utils.h"

bool dir_exists(const char *path) {
    struct stat info;

    if(stat( path, &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}
