//
// Created by sachetto on 31/05/19.
//

#include "pvd_utils.h"
#include "../xml_parser/yxml.h"
#include "../utils/file_utils.h"
#include "../single_file_libraries/stb_ds.h"

#include <stdlib.h>
#include <string.h>



static int parse_pvd_xml(yxml_t *x, yxml_ret_t r, struct vtk_files* vtk_files) {

    static char *current_file = NULL;
    static char *current_step = NULL;

    switch(r) {
        case YXML_OK:
            break;
        case YXML_ELEMSTART:
            break;
        case YXML_ELEMEND:
            break;
        case YXML_ATTRSTART:
            break;
        case YXML_ATTREND:
            if(strcmp(x->attr, "timestep") == 0) {
                arrpush(current_step, '\0');
                arrpush(vtk_files->timesteps, strtof(current_step, NULL));
                arrsetlen(current_step, 0);
            }
            else  if(strcmp(x->attr, "file") == 0) {
                arrpush(current_file, '\0');
                arrpush(vtk_files->files_list, strdup(current_file));
                arrsetlen(current_file, 0);
            }
            break;
        case YXML_PICONTENT:
        case YXML_CONTENT:
            break;
        case YXML_ATTRVAL:
            if(strcmp(x->attr, "timestep") == 0) {
                arrpush(current_step, x->data[0]);
            }
            else  if(strcmp(x->attr, "file") == 0) {
                arrpush(current_file, x->data[0]);
            }
            break;
        case YXML_PISTART:
            break;
        case YXML_PIEND:
            break;
        default:
            exit(0);
    }

    return 0;
}


struct vtk_files* list_files_from_and_timesteps_from_pvd(const char *pvd_file) {

    struct vtk_files* vtk_files = (struct vtk_files*) calloc(1, sizeof(struct vtk_files));

    if(vtk_files == NULL) {
        return NULL;
    }

    size_t size;

    char *tmp = read_entire_file_with_mmap(pvd_file, &size);

    if(!tmp) {
        return NULL;
    }

    char *source = tmp;

    //VTK XML file
    static char stack[8*1024];
    yxml_t *x = (yxml_t *) malloc(sizeof(yxml_t));
    yxml_init(x, stack, sizeof(stack));

    for (; *source; source++) {
        yxml_ret_t r = yxml_parse(x, *source);
        parse_pvd_xml(x, r, vtk_files);
    }

    free(x);
    munmap(tmp, size);

    return vtk_files;

}