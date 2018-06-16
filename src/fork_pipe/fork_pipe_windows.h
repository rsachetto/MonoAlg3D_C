//
// Created by sache on 29/05/2018.
//

#ifndef MONOALG3D_FORK_PIPE_WINDOWS_H
#define MONOALG3D_FORK_PIPE_WINDOWS_H

#include<windows.h>
#pragma comment(lib, "User32.lib")

void DisplayError(char *pszAPI, int line);
void ReadAndHandleOutput(HANDLE hOutputRead, void (*function_to_apply)(void*));
void PrepAndLaunchRedirectedChild(char * program, HANDLE hChildStdOut,
                                  HANDLE hChildStdIn,
                                  HANDLE hChildStdErr);

#endif //MONOALG3D_FORK_PIPE_WINDOWS_H
