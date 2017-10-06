//
// Created by sachetto on 03/10/17.
//

#include <sys/stat.h>
#include <malloc.h>
#include <string.h>
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


struct output_utils *new_output_utils(int print_rate, const char *output_dir) {

    struct output_utils * output_info = (struct output_utils*)malloc(sizeof(struct output_utils));
    output_info->print_rate = print_rate;
    output_info->output_dir_name = sdsnew(output_dir);

    return output_info;
}

void free_output_utils(struct output_utils* info) {
    if(info->output_dir_name) {
        sdsfree(info->output_dir_name);
    }
    free(info);
}