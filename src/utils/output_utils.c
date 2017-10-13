//
// Created by sachetto on 03/10/17.
//

#include <sys/stat.h>
#include <malloc.h>
#include <assert.h>
#include "output_utils.h"
#include "config_parser.h"

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

void configure_output_from_options(struct output_utils *output_utils,
                                              struct user_options *options) {

    assert(output_utils);
    assert(options);

    output_utils->print_rate = options->print_rate;
    output_utils->output_dir_name = sdsnew(options->out_dir_name);

}