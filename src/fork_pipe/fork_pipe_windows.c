//
// Created by sache on 29/05/2018.
//

#include "fork_pipe_windows.h"

HANDLE hStdIn = NULL; // Handle to parents std input.

//TODO: make it receive a function that will be called after reading the child output
DWORD WINAPI run_child_process_and_process_output (LPVOID lpParam)
{
    HANDLE hOutputReadTmp,hOutputRead,hOutputWrite;
    HANDLE hInputWriteTmp,hInputRead,hInputWrite;
    HANDLE hErrorWrite;
    SECURITY_ATTRIBUTES sa;

    struct ThreadData *td = (struct ThreadData*)lpParam;

    // Set up the security attributes struct.
    sa.nLength= sizeof(SECURITY_ATTRIBUTES);
    sa.lpSecurityDescriptor = NULL;
    sa.bInheritHandle = TRUE;

    // Create the child output pipe.
    if (!CreatePipe(&hOutputReadTmp,&hOutputWrite,&sa,0))
        DisplayError("CreatePipe");


    // Create a duplicate of the output write handle for the std error
    // write handle. This is necessary in case the child application
    // closes one of its std output handles.
    if (!DuplicateHandle(GetCurrentProcess(),hOutputWrite,
                         GetCurrentProcess(),&hErrorWrite,0,
                         TRUE,DUPLICATE_SAME_ACCESS))
        DisplayError("DuplicateHandle");


    // Create the child input pipe.
    if (!CreatePipe(&hInputRead,&hInputWriteTmp,&sa,0))
        DisplayError("CreatePipe");


    // Create new output read handle and the input write handles. Set
    // the Properties to FALSE. Otherwise, the child inherits the
    // properties and, as a result, non-closeable handles to the pipes
    // are created.
    if (!DuplicateHandle(GetCurrentProcess(),hOutputReadTmp,
                         GetCurrentProcess(),
                         &hOutputRead, // Address of new handle.
                         0,FALSE, // Make it uninheritable.
                         DUPLICATE_SAME_ACCESS))
        DisplayError("DupliateHandle");

    if (!DuplicateHandle(GetCurrentProcess(),hInputWriteTmp,
                         GetCurrentProcess(),
                         &hInputWrite, // Address of new handle.
                         0,FALSE, // Make it uninheritable.
                         DUPLICATE_SAME_ACCESS))
        DisplayError("DupliateHandle");


    // Close inheritable copies of the handles you do not want to be
    // inherited.
    if (!CloseHandle(hOutputReadTmp)) DisplayError("CloseHandle");
    if (!CloseHandle(hInputWriteTmp)) DisplayError("CloseHandle");


    // Get std input handle so you can close it and force the ReadFile to
    // fail when you want the input thread to exit.
    if ( (hStdIn = GetStdHandle(STD_INPUT_HANDLE)) == INVALID_HANDLE_VALUE )
        DisplayError("GetStdHandle");

    PrepAndLaunchRedirectedChild(hOutputWrite,hInputRead,hErrorWrite);


    // Close pipe handles (do not continue to modify the parent).
    // You need to make sure that no handles to the write end of the
    // output pipe are maintained in this process or else the pipe will
    // not close when the child process exits and the ReadFile will hang.
    if (!CloseHandle(hOutputWrite)) DisplayError("CloseHandle");
    if (!CloseHandle(hInputRead )) DisplayError("CloseHandle");
    if (!CloseHandle(hErrorWrite)) DisplayError("CloseHandle");

    // Read the child's output.
    ReadAndHandleOutput(hOutputRead, td->fn_pointer);
    // Redirection is complete

    // Force the read on the input to return by closing the stdin handle.
    if (!CloseHandle(hStdIn)) DisplayError("CloseHandle");

    if (!CloseHandle(hOutputRead)) DisplayError("CloseHandle");
    if (!CloseHandle(hInputWrite)) DisplayError("CloseHandle");
}


///////////////////////////////////////////////////////////////////////
// PrepAndLaunchRedirectedChild
// Sets up STARTUPINFO structure, and launches redirected child.
///////////////////////////////////////////////////////////////////////


void PrepAndLaunchRedirectedChild(HANDLE hChildStdOut,
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

    //TODO: we need to pass the child path and parameters here
    char program[] = "bin\\MonoAlg3D.exe -c example_configs\\benchmark_config_example.ini";


    //strcpy_s (program, strlen("bin\\MonoAlg3D.exe "), "bin\\MonoAlg3D.exe ");
    //strcat_s (program, strlen("-c "), "-c ");
    //strcat_s (program, strlen("example_configs/benchmark_config_example.ini"), "example_configs/benchmark_config_example.ini");

    if (!CreateProcess(NULL, program, NULL,NULL,TRUE,
                       CREATE_NEW_CONSOLE,NULL,NULL,&si,&pi))
        DisplayError("CreateProcess");

    // Close any unnecessary handles.
    if (!CloseHandle(pi.hThread)) DisplayError("CloseHandle");

    //free(program);
}


///////////////////////////////////////////////////////////////////////
// ReadAndHandleOutput
// Monitors handle for input. Exits when child exits or pipe breaks.
///////////////////////////////////////////////////////////////////////
void ReadAndHandleOutput(HANDLE hPipeRead, void (*f)(void*))
{


    char lpBuffer[256];
    DWORD nBytesRead;
    DWORD nCharsWritten;

    while(TRUE)
    {
        if (!ReadFile(hPipeRead,lpBuffer, sizeof(lpBuffer),
                      &nBytesRead,NULL) || !nBytesRead)
        {
            if (GetLastError() == ERROR_BROKEN_PIPE)
                break; // pipe done - normal exit path.
            else
                DisplayError("ReadFile"); // Something bad happened.
        }

        // Display the character read on the screen.
        //WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE),lpBuffer, nBytesRead,&nCharsWritten,NULL);

        //f will free the data passed to it
        lpBuffer[nBytesRead-1] = '\0';
        f(_strdup(lpBuffer));

    }
}

///////////////////////////////////////////////////////////////////////
// DisplayError
// Displays the error number and corresponding message.
///////////////////////////////////////////////////////////////////////
void DisplayError(char *pszAPI)
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
             "ERROR: API    = %s.\n   error code = %d.\n   message    = %s.\n",
             pszAPI, GetLastError(), (char *)lpvMessageBuffer);

    WriteConsole(GetStdHandle(STD_OUTPUT_HANDLE),szPrintBuffer,
                 lstrlen(szPrintBuffer),&nCharsWritten,NULL);

    LocalFree(lpvMessageBuffer);
    ExitProcess(GetLastError());
}