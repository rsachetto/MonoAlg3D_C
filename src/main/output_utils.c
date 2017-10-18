//
// Created by sachetto on 03/10/17.
//

#include <sys/stat.h>
#include <malloc.h>
#include <assert.h>
#include <stdlib.h>
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


struct output_utils *new_output_utils() {

    struct output_utils * output_info = (struct output_utils*)malloc(sizeof(struct output_utils));
    return output_info;
}

void free_output_utils(struct output_utils* info) {

    assert(info);

    if(info->output_dir_name) {
        sdsfree(info->output_dir_name);
    }
    free(info);
}

void create_dir_if_no_exists(const char *out_dir) {
    if (!dir_exists (out_dir)) {
        printf ("%s does not exist! Creating\n", out_dir);

        if (mkdir (out_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) {
            fprintf (stderr, "Error creating directory %s. Exiting!\n", out_dir);
            exit (EXIT_FAILURE);
        }
    }
}