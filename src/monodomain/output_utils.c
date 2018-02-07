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


bool dir_exists(const char *path) {
    struct stat info;

    if(stat( path, &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}

void create_dir_if_no_exists(const char *out_dir) {
    if (!dir_exists (out_dir)) {
        printf ("%s does not exist! Creating!\n", out_dir);


#if defined _MSC_VER
		if (_mkdir(out_dir) == -1) {
#else
		if (mkdir(out_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) {
#endif
            fprintf (stderr, "Error creating directory %s. Exiting!\n", out_dir);
            exit (EXIT_FAILURE);
        }
    }
}