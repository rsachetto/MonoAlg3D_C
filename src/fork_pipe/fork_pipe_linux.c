//
// Created by sache on 29/05/2018.
//

#include "fork_pipe_linux.h"


//TODO: make it receive a function that will be called after reading the child output
void run_child_process_and_process_output (char * program_with_path,  void (*function_to_apply)(void*)) {
    //TODO: implement fork, pipe, dup2 here



}
