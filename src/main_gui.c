#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <tchar.h>
#include <Windows.h>
#pragma warning(disable:4996)
#endif

#ifdef linux
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#include <libgen.h>

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
uiWindow *mainWindow;
uiBox *verticalBox, *horizontalButtonBox, *horizontalTextBox;
uiButton *btnConfiguration, *btnRun, *btnSave, *btnParaview, *btnCancel;
static uiProgressBar *pbar;

uiEntry *alg_3D_path_entry, *paraview_path_entry;

struct user_options *options;
bool simulation_running = false;
int modified = 0;

static char *config_file_name = NULL;
static char *global_alg3d_path = NULL;
static char *global_config_path = NULL;
static char *global_paraview_path = NULL;

#ifdef linux
pthread_t thread;
pid_t childpid;
#endif

static void append_to_output(void *data) {
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

            int progress = (int)((atof (time_string) / (options->final_time - options->dt_pde) * 100.0));

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

static void append_to_queue(void *data) {
    uiQueueMain(append_to_output, data);
}

static int on_closing_main_window(uiWindow *mw, void *data) {

    if (simulation_running) {

        int response = uiMsgBoxConfirmCancel (mw, "Confirm quit!",
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

    uiQuit ();
    return 1;
}

static int on_should_quit_main_program(void *data) {

    if (simulation_running) {

        int response = uiMsgBoxConfirmCancel (mainWindow, "Confirm quit!",
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

static int on_close_preferences_window(uiWindow *b, void *data) {
    //TODO: ask form confirmation if any preference changed
    return 1;
}

static void close_preferences_window(uiButton *b, void *data) {
    //TODO: ask form confirmation if any preference changed
    uiWindow *pref_win = uiWindow (data);
    uiControlDestroy (uiControl (pref_win));
}

static void on_open_path_clicked(uiButton *b, void *data) {
        uiEntry *entry = (uiEntry*) data;
        char *path = uiOpenFile (mainWindow);
        uiEntrySetText(entry, path);
}


static void on_open_file_clicked(uiButton *b, void *data) {

    if (config_file_name)
        uiFreeText (config_file_name);

    // TODO: implement filter for ini files
    config_file_name = uiOpenFile (mainWindow);
    if (config_file_name == NULL) {
        return;
    }

    FILE *f = fopen (config_file_name, "rb");
    fseek (f, 0, SEEK_END);
    long fsize = ftell (f);
    rewind (f);

    char *string = malloc (fsize + 1);
    fread (string, fsize, 1, f);

    if(string[fsize-1] == '\n') {
        string[fsize-1] = '\0';
    }
    else {
        string[fsize] = '\0';
    }

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
        uiMsgBoxError (mainWindow, "Invalid file!", "Error parsing ini file!");
    }
}

static void on_open_file_menu_clicked(uiMenuItem *item, uiWindow *w, void *data) {
    on_open_file_clicked(NULL, NULL);
}


static void save_simulation_config_file(uiButton *b, void *data) {

    if (!config_file_name) {
        uiMsgBoxError(mainWindow, "Error", "This should never be NULL!");
        return;
    }

    // TODO: we need to do this on windows
#ifdef linux
    uiSourceViewSaveSource (configText, config_file_name);
#endif

    uiControlDisable (uiControl (btnSave));
}

static void on_save_file_menu_clicked(uiMenuItem *item, uiWindow *w, void *data) {
    save_simulation_config_file(NULL, NULL);
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
    uiControlEnable(uiControl(btnParaview));

    free(td->program);
    free(td);

#ifdef _WIN32
    return 0;
#else
    return NULL;
#endif
}

void start_monodomain_exec () {

    struct ThreadData *td = (struct ThreadData *)malloc (sizeof (struct ThreadData));
    td->fn_pointer = append_to_queue;
    sds program = sdsnew ("");

    //TODO we need to change this acording to the user options maybe...
    program = sdscatfmt (program, "%s -c %s --draw_gl_output", global_alg3d_path, config_file_name);

#ifdef _WIN32
    HANDLE thread;
    DWORD threadId;
    td->program = _strdup (program);

    thread = CreateThread (NULL,                      // default security attributes
                           0,                         // use default stack size
                           start_program_with_thread, // thread function name
                           (LPVOID)td,                // argument to thread function
                           0,                         // use default creation flags
                           &threadId);                // returns the thread identifier

#else

    td->program = strdup (program);
    pthread_create (&thread, NULL, start_program_with_thread, (void *)td);
    pthread_detach (thread);
    sdsfree (program);

#endif

}


static void run_simulation(uiButton *b, void *data) {

    if(!global_alg3d_path) {
        uiMsgBoxError(mainWindow, "Error", "Simulator executable path is not configured!\nConfigure using edit->preferences menu!");
        return;
    }


    //TODO: config auto save on run. Make a windows version of this
#ifdef linux
    if(modified)
#else
    if(modified)
#endif
    {
        int response = uiMsgBoxConfirmCancel(mainWindow, "File not saved!", "The config file was modified. "
                                                                   "Do you want to save before running?");

        if(response == uiReturnValueCancel) {
            return;
        }
        else {
            save_simulation_config_file(NULL, NULL);
        }

    }

    uiMultilineEntrySetText (runText, "");
    uiProgressBarSetValue (pbar, 0);

    start_monodomain_exec ();
    uiControlDisable (uiControl (btnRun));
}

static void cancel_simulation(uiButton *b, void *data) {

    if (simulation_running) {

        int response = uiMsgBoxConfirmCancel (mainWindow, "Confirm cancellation!",
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
void on_simulation_config_changed(uiSourceView *e, void *data)
#else
void on_simulation_config_changed (uiMultilineEntry *e, void *data)
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

char *get_config_file_path() {

    char *gui_path, *config_path;

#ifdef linux
    size_t read_chars;
    size_t bufsize;
#endif

#ifdef _WIN32
    DWORD read_chars;
    DWORD bufsize;
#endif

    bufsize = 2048;
    gui_path = (char*) malloc(bufsize);
    config_path = (char*) malloc(strlen(".config") + bufsize);

#ifdef linux
    read_chars = readlink("/proc/self/exe", gui_path, bufsize);
#endif

#ifdef _WIN32

    read_chars = GetModuleFileName(NULL, gui_path, bufsize);
#endif

    gui_path[read_chars] = '\0';

    #ifdef linux
    char* base_dir = dirname(gui_path);
#endif

#ifdef _WIN32
    char node[ _MAX_DIR ];
    char base_dir[ _MAX_DIR ];
    char fname[ _MAX_FNAME ];
    char ext[ _MAX_EXT ];

    _splitpath( gui_path, node, base_dir, fname, ext );
    sprintf(config_path, "%s\\.config", base_dir);
#else
    sprintf(config_path, "%s/.config", base_dir);
#endif
    free(base_dir);
    return config_path;
}


void save_gui_config_file(uiButton *e, void *data) {

    char *alg_path = uiEntryText(alg_3D_path_entry);
    char *paraview_path = uiEntryText(paraview_path_entry);

    FILE *config_file = fopen(global_config_path, "w");

    if(!config_file) {
        fprintf(stderr, "Error openning %s\n", global_config_path);
        return;
    }

    fprintf(config_file, "%s\n%s\n", alg_path, paraview_path);

    if(global_alg3d_path) {
        free(global_alg3d_path);
    }

    if(global_paraview_path) {
        free(global_paraview_path);
    }

    global_alg3d_path = strdup(alg_path);
    global_paraview_path = strdup(paraview_path);

    free(alg_path);
    free(paraview_path);
    fclose(config_file);

    uiWindow *win = uiWindow (data);
    uiControlDestroy (uiControl (win));
}

static void open_configuration_window(uiMenuItem *item, uiWindow *mw, void *data) {

    uiWindow *pref_window;
    pref_window = uiNewWindow ("Preferences", 640, 150, 0);
    uiWindowSetModal(pref_window, NULL, 1);
    uiWindowSetMargined (pref_window, 1);
    uiWindowOnClosing(pref_window, on_close_preferences_window, pref_window);

    uiLabel *alg_3D_path_label = uiNewLabel("Alg 3D executable full path: ");
    alg_3D_path_entry = uiNewEntry();
    uiButton *alg_3D_path_choose = uiNewButton("Open...");
    uiButtonOnClicked(alg_3D_path_choose, on_open_path_clicked, alg_3D_path_entry);

    uiBox *alg_3D_path_button_label_box = uiNewHorizontalBox ();
    uiBoxAppend(alg_3D_path_button_label_box, uiControl(alg_3D_path_label), 0);
    uiBoxAppend(alg_3D_path_button_label_box, uiControl(alg_3D_path_entry), 1);
    uiBoxAppend(alg_3D_path_button_label_box, uiControl(alg_3D_path_choose), 0);

    uiLabel *libraries_path_label = uiNewLabel("Paraview full path: ");
    paraview_path_entry = uiNewEntry();
    uiButton *paraview_path_choose = uiNewButton("Open...");
    uiButtonOnClicked(paraview_path_choose, on_open_path_clicked, paraview_path_entry);

    uiBox *libraries_path_button_label_box = uiNewHorizontalBox ();
    uiBoxAppend(libraries_path_button_label_box, uiControl(libraries_path_label), 0);
    uiBoxAppend(libraries_path_button_label_box, uiControl(paraview_path_entry), 1);
    uiBoxAppend(libraries_path_button_label_box, uiControl(paraview_path_choose), 0);

    uiBoxSetPadded (alg_3D_path_button_label_box, 1);
    uiBoxSetPadded (libraries_path_button_label_box, 1);

    uiBox *vertical_box = uiNewVerticalBox();
    uiBoxAppend(vertical_box, uiControl(alg_3D_path_button_label_box), 0);
    uiBoxAppend(vertical_box, uiControl(libraries_path_button_label_box), 0);

    uiButton *ok = uiNewButton("Save and close");
    uiButtonOnClicked(ok, save_gui_config_file, (void*)pref_window);

    uiButton *cancel = uiNewButton("Close without saving");
    uiButtonOnClicked (cancel, close_preferences_window, (void*) pref_window);

    uiBox *ok_cancel_box = uiNewHorizontalBox();
    uiBoxAppend(ok_cancel_box, uiControl(ok), 0);
    uiBoxAppend(ok_cancel_box, uiControl(cancel), 0);
    uiBoxSetPadded (ok_cancel_box, 1);

    uiBoxSetPadded (libraries_path_button_label_box, 1);

    uiBoxAppend(vertical_box, uiControl(ok_cancel_box), 0);

    uiBoxSetPadded (vertical_box, 1);

    uiWindowSetChild (pref_window, uiControl (vertical_box));

    uiEntrySetText(alg_3D_path_entry, global_alg3d_path);
    uiEntrySetText(paraview_path_entry, global_paraview_path);

    uiControlShow (uiControl (pref_window));

}



#ifdef _WIN32
DWORD WINAPI start_paraview_with_thread(LPVOID data) {
#else
static void *start_paraview_with_thread(void *data) {
#endif
    char *program = (char*)data;

    system(program);
    sdsfree(program);

#ifdef _WIN32
    return 1;
#else
    return NULL;
#endif

}

void open_paraview_last_simulation(uiButton *b, void *data) {

    struct ThreadData *td = (struct ThreadData *)malloc (sizeof (struct ThreadData));
    td->fn_pointer = append_to_queue;
    sds program = sdsnew ("");

    if(options->save_mesh_config->out_dir_name == NULL) {
        uiMsgBoxError(mainWindow, "Error", "The last simulation was not saved!");
        return;
    }


#ifdef _WIN32
    HANDLE thread;
    DWORD threadId;

    program = sdscatfmt (program, "%s --data=%s\\V_t_..vtk", global_paraview_path, options->save_mesh_config->out_dir_name);
    td->program = _strdup (program);

    thread = CreateThread (NULL,                      // default security attributes
                           0,                         // use default stack size
                           start_paraview_with_thread, // thread function name
                           (LPVOID)program,                // argument to thread function
                           0,                         // use default creation flags
                           &threadId);                // returns the thread identifier

#else
    program = sdscatfmt (program, "%s --data=%s/V_t_..vtk", global_paraview_path, options->save_mesh_config->out_dir_name);
    td->program = strdup (program);
    pthread_create (&thread, NULL, start_paraview_with_thread, (void *)program);
    pthread_detach (thread);
#endif

}

void set_global_paths() {
    global_config_path = get_config_file_path();
    FILE *config_file = fopen(global_config_path, "r");

    if(config_file) {

        char *alg_path = (char*) malloc(2048);
        char *paraview_path = (char*) malloc(2048);

        fscanf(config_file, "%s\n%s\n", alg_path, paraview_path);

        char *pos;
        if ((pos=strchr(alg_path, '\n')) != NULL)
            *pos = '\0';

        if ((pos=strchr(paraview_path, '\n')) != NULL)
            *pos = '\0';

        global_alg3d_path = strdup(alg_path);
        global_paraview_path = strdup(paraview_path);

        free(alg_path);
        free(paraview_path);
        return;
    }

}

int main (int argc, char **argv) {

    options = new_user_options ();
    set_global_paths();

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
    uiMenuItemOnClicked(item, on_open_file_menu_clicked, NULL);

    item = uiMenuAppendItem (menu, "Save");
    uiMenuItemOnClicked(item, on_save_file_menu_clicked, NULL);

    item = uiMenuAppendQuitItem (menu);

    menu = uiNewMenu ("Edit");
    item = uiMenuAppendPreferencesItem (menu);
    uiMenuItemOnClicked(item, open_configuration_window, NULL);

    mainWindow = uiNewWindow ("MonoAlg3D GUI", 1024, 768, 1);
    uiWindowSetMargined (mainWindow, 1);

    horizontalButtonBox = uiNewHorizontalBox ();
    uiBoxSetPadded (horizontalButtonBox, 1);
    uiWindowSetChild (mainWindow, uiControl (horizontalButtonBox));

    horizontalTextBox = uiNewHorizontalBox ();
    uiBoxSetPadded (horizontalTextBox, 1);
    uiWindowSetChild (mainWindow, uiControl (horizontalTextBox));

    verticalBox = uiNewVerticalBox ();
    uiBoxSetPadded (verticalBox, 1);
    uiWindowSetChild (mainWindow, uiControl (verticalBox));

#ifdef linux
    configText = uiNewSourceView ();
    uiSourceViewSetReadOnly (configText, true);
    uiSourceViewOnChanged(configText, on_simulation_config_changed, NULL);
#endif

#ifdef _WIN32
    configText = uiNewMultilineEntry ();
    uiMultilineEntrySetReadOnly (configText, 1);
    uiMultilineEntryOnChanged(configText, on_simulation_config_changed, NULL);
#endif

    runText = uiNewMultilineEntry ();
    uiMultilineEntrySetReadOnly (runText, 1);

    btnSave = uiNewButton ("Save");
    uiControlDisable (uiControl (btnSave));
    uiButtonOnClicked(btnSave, save_simulation_config_file, NULL);

    btnConfiguration = uiNewButton ("Open configuration");
    uiButtonOnClicked(btnConfiguration, on_open_file_clicked, NULL);

    btnRun = uiNewButton ("Run simulation");
    uiButtonOnClicked(btnRun, run_simulation, NULL);
    uiControlDisable (uiControl (btnRun));

    btnCancel = uiNewButton ("Cancel Simulation");
    uiButtonOnClicked(btnCancel, cancel_simulation, NULL);
    uiControlDisable (uiControl (btnCancel));

    btnParaview = uiNewButton ("Visualize results");
    uiButtonOnClicked(btnParaview, open_paraview_last_simulation, NULL);
    uiControlDisable (uiControl (btnParaview));

    pbar = uiNewProgressBar ();

    uiBoxAppend (horizontalTextBox, uiControl (configText), 1);
    uiBoxAppend (horizontalTextBox, uiControl (runText), 1);

    uiBoxAppend (horizontalButtonBox, uiControl (btnConfiguration), 0);
    uiBoxAppend (horizontalButtonBox, uiControl (btnSave), 0);
    uiBoxAppend (horizontalButtonBox, uiControl (btnRun), 0);
    uiBoxAppend (horizontalButtonBox, uiControl (btnParaview), 0);
    uiBoxAppend (horizontalButtonBox, uiControl (pbar), 1);
    uiBoxAppend (horizontalButtonBox, uiControl (btnCancel), 0);

    uiBoxAppend (verticalBox, uiControl (horizontalTextBox), 1);
    uiBoxAppend (verticalBox, uiControl (horizontalButtonBox), 0);

    uiWindowOnClosing(mainWindow, on_closing_main_window, NULL);
    uiOnShouldQuit(on_should_quit_main_program, mainWindow);
    uiControlShow (uiControl (mainWindow));
    uiMain ();
    uiUninit();

    return 0;
}
