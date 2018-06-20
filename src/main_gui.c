#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <tchar.h>
#include <windows.h>
#endif

#ifdef linux
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#endif

#include "config/config_parser.h"
#include "fork_pipe/fork_pipe.h"
#include "ini_parser/ini.h"
#include "libui/ui.h"
#include "string/sds.h"

#ifdef linux
uiSourceView *configText;
#endif

#ifdef _WIN32
uiMultilineEntry *configText;
#endif

uiMultilineEntry *runText;
uiWindow *w;
uiBox *verticalBox, *horizontalButtonBox, *horizontalTextBox;
uiButton *btnConfiguration, *btnRun, *btnSave, *btnCancel;
static uiProgressBar *pbar;

struct user_options *options;
bool simulation_running = false;
int modified;

char *config_file_name = NULL;

#ifdef linux
pthread_t thread;
pid_t childpid;
#endif

static void appendToOutput (void *data) {
    char *s = (char *)data;
    uiMultilineEntryAppend (runText, s);

    char *sub = strstr (s, "t = ");

    if (sub != NULL) {

        char *sub2 = strstr (sub, ",");

        if (sub2 != NULL) {

            char *e = strchr (sub, ',');
            int index = (int)(e - sub);

            char *time_string = (char *)malloc (index - 3);

            strncpy (time_string, sub + 3, index - 3);
            time_string[index - 3] = '\0';

            int progress = (int)((atof (time_string) / (options->final_time - options->dt_edp) * 100.0));

            uiProgressBarSetValue (pbar, progress);
        }
    }

    char *sub3 = strstr (s, "CG Total Iterations: ");

    if (sub3) {
        uiProgressBarSetValue (pbar, 100);
        uiControlEnable (uiControl (btnRun));
    }

    free (s);
}

static void appendToQueue (void *data) {
    uiQueueMain (appendToOutput, data);
}

static int onClosing (uiWindow *w, void *data) {

    if (simulation_running) {

        int response = uiMsgBoxConfirmCancel (w, "Confirm quit!",
                                              "Are you sure you want to quit and cancel the current simulation?");

        if (response == uiReturnValueCancel) {
            printf ("NO CLOSE\n");
            return 0;
        } else {
#ifdef _WIN32
            // TODO: implement win32 call
#else
            if (childpid != -1) {
                kill (childpid, SIGKILL);
            }
#endif
        }
    }

    uiQuit ();
    return 1;
}

static int onShouldQuit (void *data) {

    if (simulation_running) {

        int response = uiMsgBoxConfirmCancel (w, "Confirm quit!",
                                              "Are you sure you want to quit and cancel the current simulation?");

        if (response == uiReturnValueCancel) {
            return 0;
        } else {
#ifdef _WIN32
            // TODO: implement win32 call
#else
            if (childpid != -1) {
                kill (childpid, SIGKILL);
            }
#endif
        }
    }

    uiWindow *mainwin = uiWindow (data);
    uiControlDestroy (uiControl (mainwin));
    return 1;
}

static void onOpenFileClicked (uiButton *b, void *data) {

    if (config_file_name)
        uiFreeText (config_file_name);

    // TODO: implement filter for ini files
    config_file_name = uiOpenFile (w);
    if (config_file_name == NULL) {
        return;
    }

    FILE *f = fopen (config_file_name, "rb");
    fseek (f, 0, SEEK_END);
    long fsize = ftell (f);
    rewind (f);

    char *string = malloc (fsize + 1);
    fread (string, fsize, 1, f);
#ifdef linux
    string[fsize-1] = '\0';
#else
    string[fsize] = '\0';
#endif
    fclose (f);

    if (ini_parse (config_file_name, parse_config_file, options) < 0) {
        fprintf (stderr, "Error: Can't load the config file %s\n", config_file_name);
    }

    // TODO: we should check for more missing parameters
    if (options->main_found) {

        uiControlEnable (uiControl (btnRun));

#ifdef linux
        uiSourceViewSetText (configText, string, "text/x-ini-file");
        uiSourceViewSetReadOnly (configText, false);
#endif

#ifdef _WIN32
        uiMultilineEntryAppend (configText, string);
        uiMultilineEntrySetReadOnly (configText, false);
#endif
    } else {
        uiMsgBoxError (w, "Invalid file!", "Error parsing ini file!");
    }
}

static void saveConfigFile () {

    if (!config_file_name) {
        uiMsgBoxError(w, "Error", "This should never be NULL!");
    }

    // TODO: we need to do this on windows
#ifdef linux
    uiSourceViewSaveSource (configText, config_file_name);
#endif

    uiControlDisable (uiControl (btnSave));
}

static void saveConfigFileCB (uiButton *b, void *data) {
    saveConfigFile ();
}

struct ThreadData {
    char *program;
    void (*fn_pointer) (void *);
};

#ifdef _WIN32
DWORD WINAPI start_program_with_thread (LPVOID thread_param) {
#else
void *start_program_with_thread (void *thread_param) {
#endif

    struct ThreadData *td = (struct ThreadData *)thread_param;
    simulation_running = true;
    uiControlEnable (uiControl (btnCancel));
#ifdef _WIN32
    run_child_process_and_process_output (td->program, td->fn_pointer);
#else
    run_child_process_and_process_output (td->program, td->fn_pointer, &childpid);
#endif

    simulation_running = false;
    uiControlDisable (uiControl (btnCancel));


    free(td->program);
    free(td);


    return NULL;
}

void start_monodomain_exec () {

    // LEAK: maybe free it on the run_child_process_and_process_output function
    struct ThreadData *td = (struct ThreadData *)malloc (sizeof (struct ThreadData));
    td->fn_pointer = appendToQueue;

#ifdef _WIN32
    HANDLE thread;
    DWORD threadId;

    // TODO: we need to pass the Monoalg path and parameters here
    char program[] = "bin\\MonoAlg3D.exe -c example_configs\\benchmark_config_example.ini";
    td->program = _strdup (program);

    thread = CreateThread (NULL,                      // default security attributes
                           0,                         // use default stack size
                           start_program_with_thread, // thread function name
                           (LPVOID)td,                // argument to thread function
                           0,                         // use default creation flags
                           &threadId);                // returns the thread identifier

#else
    // TODO: we need to pass the Monoalg path and parameters here

    sds program = sdsnew ("");
    program = sdscatfmt (program, "%bin/MonoAlg3D -c %s", config_file_name);
    td->program = strdup (program);
    pthread_create (&thread, NULL, start_program_with_thread, (void *)td);
    pthread_detach (thread);
    sdsfree (program);

#endif
}

// TODO: the working directory and the executable need to be read from a configuration file
static void runSimulation (uiButton *b, void *data) {

    //TODO: config auto save on run. Make a windows version of this
#ifdef linux
    if(uiSourceViewGetModified(configText)) {
        int response = uiMsgBoxConfirmCancel(w, "File not saved!", "The config file was modified. "
                                                                   "Do you want to save before running?");

        if(response == uiReturnValueCancel) {
            return;
        }
        else {
           printf("SAVING FILE\n");
           saveConfigFile();
        }

    }
#else
    if(uiMultilineEntryGetModified(configText)) {
        int response = uiMsgBoxConfirmCancel(w, "File not saved!", "The config file was modified. "
                                                                   "Do you want to save before running?");

        if(response == uiReturnValueCancel) {
            return;
        }
        else {
            printf("SAVING FILE\n");
            saveConfigFile();
        }

    }
#endif
    uiMultilineEntrySetText (runText, "");
    uiProgressBarSetValue (pbar, 0);

    start_monodomain_exec ();
    uiControlDisable (uiControl (btnRun));
}

static void cancelSimulation (uiButton *b, void *data) {

    if (simulation_running) {

        int response = uiMsgBoxConfirmCancel (w, "Confirm cancellation!",
                                              "Are you sure you want to cancel the current simulation?");

        if (response == uiReturnValueCancel) {
            return;
        } else {
#ifdef _WIN32
            // TODO: implement win32 call
#else
            if (childpid != -1) {
                kill (childpid, SIGKILL);
            }
#endif

            uiControlDisable (uiControl (btnCancel));
            uiProgressBarSetValue (pbar, 0);
            uiControlEnable (uiControl (btnRun));
            simulation_running = false;
        }
    }
}

#ifdef linux
void onConfigChanged (uiSourceView *e, void *data)
#else
void onConfigChanged (uiMultilineEntry *e, void *data)
#endif
{

    // TODO: make a windows version of this function

#ifdef linux
    modified = uiSourceViewGetModified (e);

    if (modified) {
        uiControlEnable (uiControl (btnSave));
    } else {
        uiControlDisable (uiControl (btnSave));
    }
#else

    modified = 1;
    uiControlEnable (uiControl (btnSave));

#endif
}

int main () {

    options = new_user_options ();

    uiMenu *menu;
    uiMenuItem *item;

    uiInitOptions o;

    memset(&o, 0, sizeof (uiInitOptions));
    const char *err = uiInit(&o);
    if (err != NULL) {
        fprintf(stderr, "error initializing ui: %s\n", err);
        uiFreeInitError(err);
        return 1;
    }


    menu = uiNewMenu ("File");
    item = uiMenuAppendItem (menu, "Open");
    // uiMenuItemOnClicked(item, openClicked, NULL);

    item = uiMenuAppendItem (menu, "Save");
    // uiMenuItemOnClicked(item, saveClicked, NULL);

    item = uiMenuAppendQuitItem (menu);

    menu = uiNewMenu ("Edit");
    item = uiMenuAppendPreferencesItem (menu);
    // uiMenuItemOnClicked(item, saveClicked, NULL);

    w = uiNewWindow ("MonoAlg3D GUI", 1024, 768, 1);
    uiWindowSetMargined (w, 1);

    horizontalButtonBox = uiNewHorizontalBox ();
    uiBoxSetPadded (horizontalButtonBox, 1);
    uiWindowSetChild (w, uiControl (horizontalButtonBox));

    horizontalTextBox = uiNewHorizontalBox ();
    uiBoxSetPadded (horizontalTextBox, 1);
    uiWindowSetChild (w, uiControl (horizontalTextBox));

    verticalBox = uiNewVerticalBox ();
    uiBoxSetPadded (verticalBox, 1);
    uiWindowSetChild (w, uiControl (verticalBox));

#ifdef linux
    configText = uiNewSourceView ();
    uiSourceViewSetReadOnly (configText, true);
    uiSourceViewOnChanged (configText, onConfigChanged, NULL);
#endif

#ifdef _WIN32
    configText = uiNewMultilineEntry ();
    uiMultilineEntrySetReadOnly (configText, 1);
    uiMultilineEntryOnChanged(configText, onConfigChanged, NULL);
#endif

    runText = uiNewMultilineEntry ();
    uiMultilineEntrySetReadOnly (runText, 1);

    btnSave = uiNewButton ("Save");
    uiControlDisable (uiControl (btnSave));
    uiButtonOnClicked (btnSave, saveConfigFileCB, NULL);

    btnConfiguration = uiNewButton ("Open configuration");
    uiButtonOnClicked (btnConfiguration, onOpenFileClicked, NULL);

    btnRun = uiNewButton ("Run simulation");
    uiButtonOnClicked (btnRun, runSimulation, NULL);
    uiControlDisable (uiControl (btnRun));

    btnCancel = uiNewButton ("Cancel Simulation");
    uiButtonOnClicked (btnCancel, cancelSimulation, NULL);
    uiControlDisable (uiControl (btnCancel));

    pbar = uiNewProgressBar ();

    uiBoxAppend (horizontalTextBox, uiControl (configText), 1);
    uiBoxAppend (horizontalTextBox, uiControl (runText), 1);

    uiBoxAppend (horizontalButtonBox, uiControl (btnConfiguration), 0);
    uiBoxAppend (horizontalButtonBox, uiControl (btnSave), 0);
    uiBoxAppend (horizontalButtonBox, uiControl (btnRun), 0);
    uiBoxAppend (horizontalButtonBox, uiControl (pbar), 1);
    uiBoxAppend (horizontalButtonBox, uiControl (btnCancel), 0);

    uiBoxAppend (verticalBox, uiControl (horizontalTextBox), 1);
    uiBoxAppend (verticalBox, uiControl (horizontalButtonBox), 0);

    uiWindowOnClosing (w, onClosing, NULL);
    uiOnShouldQuit (onShouldQuit, w);
    uiControlShow (uiControl (w));
    uiMain ();
    uiUninit();

    return 0;
}
