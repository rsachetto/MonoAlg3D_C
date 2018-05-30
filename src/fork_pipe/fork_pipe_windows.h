//
// Created by sache on 29/05/2018.
//

#ifndef MONOALG3D_FORK_PIPE_WINDOWS_H
#define MONOALG3D_FORK_PIPE_WINDOWS_H

#include<windows.h>
#pragma comment(lib, "User32.lib")

DWORD WINAPI run_child_process_and_process_output (LPVOID params);
void DisplayError(char *pszAPI);
void ReadAndHandleOutput(HANDLE hOutputRead, void (*f)(void*));
void PrepAndLaunchRedirectedChild(HANDLE hChildStdOut,
                                  HANDLE hChildStdIn,
                                  HANDLE hChildStdErr);


struct ThreadData {
    char *program;
    void (*fn_pointer)(void*);
};

#endif //MONOALG3D_FORK_PIPE_WINDOWS_H
