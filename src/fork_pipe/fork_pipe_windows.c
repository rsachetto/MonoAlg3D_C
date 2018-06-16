//
// Created by sache on 29/05/2018.
//

#include <assert.h>
#include <stdio.h>
#include "fork_pipe_windows.h"

void run_child_process_and_process_output (char * program_with_path,  void (*function_to_apply)(void*))
{
    HANDLE hStdIn = NULL; // Handle to parents std input.
    HANDLE hOutputReadTmp, hOutputRead, hOutputWrite;
    HANDLE hInputWriteTmp,hInputRead,hInputWrite;
    HANDLE hErrorWrite;
    SECURITY_ATTRIBUTES sa;

    // Set up the security attributes struct.
    sa.nLength= sizeof(SECURITY_ATTRIBUTES);
    sa.lpSecurityDescriptor = NULL;
    sa.bInheritHandle = TRUE;

    // Create the child output pipe.
    if (!CreatePipe(&hOutputReadTmp,&hOutputWrite,&sa,0))
        DisplayError("CreatePipe", __LINE__);

    // Create a duplicate of the output write handle for the std error
    // write handle. This is necessary in case the child application
    // closes one of its std output handles.
    if (!DuplicateHandle(GetCurrentProcess(),hOutputWrite,
                         GetCurrentProcess(),&hErrorWrite,0,
                         TRUE,DUPLICATE_SAME_ACCESS))
        DisplayError("DuplicateHandle", __LINE__);


    // Create the child input pipe.
    if (!CreatePipe(&hInputRead,&hInputWriteTmp,&sa,0))
        DisplayError("CreatePipe", __LINE__);


    // Create new output read handle and the input write handles. Set
    // the Properties to FALSE. Otherwise, the child inherits the
    // properties and, as a result, non-closeable handles to the pipes
    // are created.
    if (!DuplicateHandle(GetCurrentProcess(),hOutputReadTmp,
                         GetCurrentProcess(),
                         &hOutputRead, // Address of new handle.
                         0,FALSE, // Make it uninheritable.
                         DUPLICATE_SAME_ACCESS))
        DisplayError("DuplicateHandle", __LINE__);

    if (!DuplicateHandle(GetCurrentProcess(),hInputWriteTmp,
                         GetCurrentProcess(),
                         &hInputWrite, // Address of new handle.
                         0,FALSE, // Make it uninheritable.
                         DUPLICATE_SAME_ACCESS))
        DisplayError("DuplicateHandle", __LINE__);


    // Close inheritable copies of the handles you do not want to be
    // inherited.
    if (!CloseHandle(hOutputReadTmp)) DisplayError("CloseHandle", __LINE__);
    if (!CloseHandle(hInputWriteTmp)) DisplayError("CloseHandle", __LINE__);

    // Get std input handle so you can close it and force the ReadFile to
    // fail when you want the input thread to exit.
    if ((hStdIn = GetStdHandle(STD_INPUT_HANDLE)) == INVALID_HANDLE_VALUE) {
        DisplayError("GetStdHandle", __LINE__);
    }

    PrepAndLaunchRedirectedChild(program_with_path, hOutputWrite,hInputRead,hErrorWrite);

    // Close pipe handles (do not continue to modify the parent).
    // You need to make sure that no handles to the write end of the
    // output pipe are maintained in this process or else the pipe will
    // not close when the child process exits and the ReadFile will hang.
    if (!CloseHandle(hOutputWrite)) DisplayError("CloseHandle", __LINE__);
    if (!CloseHandle(hInputRead )) DisplayError("CloseHandle", __LINE__);
    if (!CloseHandle(hErrorWrite)) DisplayError("CloseHandle", __LINE__);

    // Read the child's output.
    ReadAndHandleOutput(hOutputRead, function_to_apply);
    // Redirection is complete

    // Force the read on the input to return by closing the stdin handle.
    //TODO: BUG, if we run this fork multiple times, close handle returns a error!
    // error code = 6.
    //if (!CloseHandle(hStdIn)) DisplayError("CloseHandle", __LINE__);

    if (!CloseHandle(hOutputRead)) DisplayError("CloseHandle",__LINE__ );
    if (!CloseHandle(hInputWrite)) DisplayError("CloseHandle", __LINE__);
}


///////////////////////////////////////////////////////////////////////
// PrepAndLaunchRedirectedChild
// Sets up STARTUPINFO structure, and launches redirected child.
///////////////////////////////////////////////////////////////////////


void PrepAndLaunchRedirectedChild(char *program, HANDLE hChildStdOut,
                                  HANDLE hChildStdIn,
                                  HANDLE hChildStdErr)
{
    PROCESS_INFORMATION pi;
    STARTUPINFO si;

    // Set up the start up info struct.
    ZeroMemory(&si,sizeof(STARTUPINFO));
    si.cb = sizeof(STARTUPINFO);
    si.dwFlags = STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW;
    si.hStdOutput = hChildStdOut;
    si.hStdInput  = hChildStdIn;
    si.hStdError  = hChildStdErr;
    // Use this if you want to hide the child:
    si.wShowWindow = SW_HIDE;
    // Note that dwFlags must include STARTF_USESHOWWINDOW if you want to
    // use the wShowWindow flags.


    // Launch the process that you want to redirect (in this case,
    // Child.exe). Make sure Child.exe is in the same directory as
    // redirect.c launch redirect from a command line to prevent location
    // confusion.
    if (!CreateProcess(NULL, program, NULL,NULL,TRUE,
                       CREATE_NEW_CONSOLE,NULL,NULL,&si,&pi))
        DisplayError("CreateProcess", __LINE__);

    // Close any unnecessary handles.
    if (!CloseHandle(pi.hThread)) DisplayError("CloseHandle", __LINE__);

}


///////////////////////////////////////////////////////////////////////
// ReadAndHandleOutput
// Monitors handle for input. Exits when child exits or pipe breaks.
///////////////////////////////////////////////////////////////////////
void ReadAndHandleOutput(HANDLE hPipeRead, void (*function_to_apply)(void*))
{
    char lpBuffer[1024];
    char ch;
    DWORD nBytesRead;
    DWORD nCharsWritten;

    int char_count = 0;

    while(TRUE)
    {
        if (!ReadFile(hPipeRead, &ch, 1,
                      &nBytesRead,NULL) || !nBytesRead)
        {
            if (GetLastError() == ERROR_BROKEN_PIPE)
                break; // pipe done - normal exit path.
            else
                DisplayError("ReadFile",__LINE__); // Something bad happened.
        }

        if(ch != '\n') {
            lpBuffer[char_count] = ch;
            char_count++;
        } else {
            assert(char_count < 1024);
            lpBuffer[char_count] = '\n';
            lpBuffer[char_count+1] = '\0';
            char_count = 0;
            function_to_apply(_strdup(lpBuffer));
        }

    }
}

///////////////////////////////////////////////////////////////////////
// DisplayError
// Displays the error number and corresponding message.
///////////////////////////////////////////////////////////////////////
void DisplayError(char *pszAPI, int line)
{
    LPVOID lpvMessageBuffer;
    CHAR szPrintBuffer[512];
    DWORD nCharsWritten;

    FormatMessage(
            FORMAT_MESSAGE_ALLOCATE_BUFFER|FORMAT_MESSAGE_FROM_SYSTEM,
            NULL, GetLastError(),
            MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
            (LPTSTR)&lpvMessageBuffer, 0, NULL);

    wsprintf(szPrintBuffer,
             "LINE:%d, ERROR: API    = %s.\n   error code = %d.\n   message    = %s.\n",
             line, pszAPI, GetLastError(), (char *)lpvMessageBuffer);

    WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE),szPrintBuffer,
                 lstrlen(szPrintBuffer),&nCharsWritten,NULL);

    LocalFree(lpvMessageBuffer);
    //ExitProcess(GetLastError());
}