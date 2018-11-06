//
// Created by sache on 29/05/2018.
//

#include <string.h>
#include <assert.h>
#include <fcntl.h>
#include "fork_pipe_linux.h"

void run_child_process_and_process_output (char * program_with_path,  void (*function_to_apply)(void*), pid_t *childpid) {

    int cp[2]; //child to parent pipe
    *childpid = -1;

    pipe(cp);

    char *token;


    //TODO: we need to pass more args here
    char **prog = (char**) malloc(5*sizeof(char*));

    //get the first token
    token = strtok(program_with_path, " ");

    int i = 0;
    while( token != NULL ) {

        prog[i] = strdup(token);

        token = strtok(NULL, " ");

        i++;
        //assert(i <= 4);
    }

    prog[4] = NULL;

    char ch;

    int MAX_EXPECTED_LINE = 1024;
    char line[MAX_EXPECTED_LINE];
    int char_count = 0;

    switch (*childpid = fork() ) {
        case -1: {
            perror("fork");
            exit(EXIT_FAILURE);
        }

        case 0: {

            close(cp[0]);
            dup2(cp[1], 1);

            execvp(prog[0], prog);
            perror("Error on execvp");
            exit(EXIT_FAILURE);
        }

        default: {
            close(cp[1]);
            ssize_t r;
            r = read(cp[0], &ch, 1);

            while(r == 1) {

                if(ch != '\n') {
                    line[char_count] = ch;
                    char_count++;
                } else {
                    assert(char_count < MAX_EXPECTED_LINE);
                    line[char_count] = '\n';
                    line[char_count+1] = '\0';
                    char_count = 0;
                    function_to_apply(strdup(line));
                }

                r = read(cp[0], &ch, 1);

            }

            *childpid = -1;

        }
    }



}
