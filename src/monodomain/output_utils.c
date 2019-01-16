//
// Created by sachetto on 03/10/17.
//

#if defined _MSC_VER
#include <direct.h>
#endif

#include "output_utils.h"
#include <sys/stat.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>


bool dir_exists(const char *path) {
    struct stat info;

    if(stat( path, &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}

void create_dir(const char *out_dir) {

    //TODO: check for windows dir separators
    int dirs_count;

    sds *all_dirs = sdssplit(out_dir, "/", &dirs_count);
    sds new_dir = sdsempty();

    for(int d = 0; d < dirs_count; d++) {

        new_dir = sdscat(new_dir, all_dirs[d]);
        new_dir = sdscat(new_dir, "/");

        if (!dir_exists (new_dir)) {

            printf ("%s does not exist! Creating!\n", new_dir);
            #if defined _MSC_VER
                if (_mkdir(out_dir) == -1)
            #else
                if (mkdir(new_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
            #endif
            {

                fprintf (stderr, "Error creating directory %s Exiting!\n", new_dir);
                exit (EXIT_FAILURE);
            }
        }
    }
    
}
