#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <stdbool.h>
#include "ini_parser/ini.h"
#include "config/config_parser.h"

#ifdef _WIN32
#include <windows.h>
#include "fork_pipe/fork_pipe_windows.h" //TODO: put this in one header
#include <tchar.h>
#endif

#include "libui/ui.h"


#ifdef linux
#include <unistd.h>
#endif

//#include <glib.h> We are trying to substitute glib for now

#ifdef linux
uiSourceView *configText;
#endif

#ifdef _WIN32
uiMultilineEntry *configText;
#endif

uiMultilineEntry *runText;
uiWindow *w;
uiBox *verticalBox, *horizontalButtonBox, *horizontalTextBox;
uiButton *btnConfiguration, *btnRun, *btnSave, *btnSaveAs;
static uiProgressBar *pbar;

struct user_options *options;
bool child_running = false;

char *config_file_name = NULL;
char *output_last_sim = NULL;

/*
gint child_stdout, child_stderr;

GIOChannel *channel_stderr, *channel_stdout;
GPid child_pid;
GThread *thread;
*/

void appendToOutput(void *data) {
    char *s = (char*) data;
    uiMultilineEntryAppend(runText, s);

    char *sub = strstr(s, "t = ");

    if(sub != NULL) {

        printf(sub);
        char *e = strchr(sub, ',');
        int index = (int)(e - sub);

        char *time_string = (char*) malloc(index-3);

        strncpy(time_string, sub+3, index-3);
        time_string[index-3] = '\0';

        int progress = (int) ((atof(time_string) / (options->final_time-options->dt_edp)*100.0));

        uiProgressBarSetValue(pbar, progress);
    }

    free(s);
}

void appendToQueue(void *data) {
    uiQueueMain(appendToOutput, data);
}


int onClosing(uiWindow *w, void *data)
{

/*
    if(child_running) {
        g_io_channel_shutdown(channel_stderr, TRUE, NULL);
        g_io_channel_shutdown(channel_stdout, TRUE, NULL);
        g_spawn_close_pid(child_pid);
    }
*/

    uiQuit();
    return 1;
}

static void onOpenFileClicked(uiButton *b, void *data)
{

    if(config_file_name) uiFreeText(config_file_name);

    //TODO: implement filter for ini files
    config_file_name = uiOpenFile(w);
    if (config_file_name == NULL) {
        return;
    }

    FILE *f = fopen(config_file_name, "rb");
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    rewind(f);

    char *string = malloc(fsize + 1);
    fread(string, fsize, 1, f);
    fclose(f);

    if (ini_parse (config_file_name, parse_config_file, options) < 0) {
        fprintf (stderr, "Error: Can't load the config file %s\n", config_file_name);
    }

    //TODO: we should check for more missing parameters
    if(options->main_found) {
		
		uiControlEnable(uiControl(btnRun));
		
		#ifdef linux
        uiSourceViewSetText(configText, string, "text/x-ini-file");
		#endif        
		
		#ifdef linux
        uiSourceViewSetReadOnly(configText, false);
		#endif
		
		#ifdef _WIN32
        uiMultilineEntryAppend(configText, string);
		uiMultilineEntrySetReadOnly(configText, false);
		#endif
    }
    else {
        uiMsgBoxError(w, "Invalid file!", "Error parsing ini file!");
    }
}

static void saveConfigFile()
{
    if(!config_file_name) uiMsgBoxError(w, "Error", "This should never be NULL!");

	
	//TODO: we need to do this on windows
	
	#ifdef linux
    uiSourceViewSaveSource(configText, config_file_name);
	#endif
    uiControlDisable(uiControl(btnSave));
}

static void saveConfigFileCB(uiButton *b, void *data)
{
   saveConfigFile();
}

/*
static gboolean getOutput (GIOChannel *channel, GIOCondition cond, gpointer data) {
    gchar *str_return;
    gsize length;
    gsize terminator_pos;
    GError *error = NULL;


    if (cond & G_IO_HUP) return FALSE;

    if (g_io_channel_read_line (channel, &str_return, &length, &terminator_pos, &error) == G_IO_STATUS_ERROR) {
        g_warning ("Something went wrong");
    }
    if (error != NULL) {
        g_warning (error->message);
    }

    uiMultilineEntryAppend(runText, str_return);

    char *sub = strstr(str_return, "t = ");

    if(sub != NULL) {

        printf(sub);
        char *e = strchr(sub, ',');
        int index = (int)(e - sub);

        char *time_string = (char*) malloc(index-3);

        strncpy(time_string, sub+3, index-3);
        time_string[index-3] = '\0';

        int progress = (int) ((atof(time_string) / (options->final_time-options->dt_edp)*100.0));

        uiProgressBarSetValue(pbar, progress);
    }

    uiHandlePendingEvents();

    g_free(str_return);

    return TRUE;
}

static gboolean getError (GIOChannel *channel, GIOCondition cond, gpointer data) {
    gchar *str_return;
    gsize length;
    gsize terminator_pos;
    GError *error = NULL;

    if (cond & G_IO_HUP) return FALSE;

    if (g_io_channel_read_line (channel, &str_return, &length, &terminator_pos, &error) == G_IO_STATUS_ERROR) {
        g_warning ("Something went wrong");
    }
    if (error != NULL) {
        g_warning (error->message);
    }


    uiMultilineEntryAppend(runText, str_return);
    uiHandlePendingEvents();

    g_free(str_return);

    return TRUE;
}


static void child_watch_cb (GPid pid, gint status, gpointer user_data) {
    //g_message ("Child %" G_PID_FORMAT " exited %s", pid, g_spawn_check_exit_status (status, NULL) ? "normally" : "abnormally");

    g_io_channel_shutdown(channel_stderr,TRUE,NULL);
    g_io_channel_shutdown(channel_stdout,TRUE,NULL);
    g_spawn_close_pid (pid);
    child_running = false;
    uiControlEnable(uiControl(btnRun));
}
*/

void start_monodomain_exec() {

    char **program = (char**) malloc(4*sizeof(char*));

    program[0] = strdup("bin/MonoAlg3D");
    program[1] = strdup("-c");
    program[2] = strdup(config_file_name);
    printf ("%s\n", config_file_name);

    program[3] = NULL;


    struct ThreadData *td = (struct ThreadData *)malloc(sizeof(struct ThreadData));

    td->program = NULL;
    td->fn_pointer = appendToQueue;

#ifdef _WIN32
    HANDLE thread;
    DWORD  threadId;
    thread = CreateThread(
            NULL,                   // default security attributes
            0,                      // use default stack size
            run_child_process_and_process_output,       // thread function name
            (LPVOID) td,          // argument to thread function
            0,                      // use default creation flags
            &threadId);   // returns the thread identifier
#endif
	/*
    g_autoptr(GError) error = NULL;

    // Spawn child process.
    g_spawn_async_with_pipes ("/home/sachetto/Projects/MonoAlg3D_C/", program, NULL, G_SPAWN_DO_NOT_REAP_CHILD, NULL,
                              NULL, &child_pid, NULL, &child_stdout,
                              &child_stderr, &error);
    if (error != NULL)
    {
        g_error ("Spawning child failed: %s", error->message);
    }

    child_running = true;
//
    g_child_watch_add (child_pid, child_watch_cb, NULL);

    channel_stderr = g_io_channel_unix_new(child_stderr);
    g_io_add_watch (channel_stderr, G_IO_IN | G_IO_HUP, (GIOFunc) getError, NULL);

    channel_stdout = g_io_channel_unix_new(child_stdout);
    g_io_add_watch (channel_stdout, G_IO_IN | G_IO_HUP, (GIOFunc) getOutput, NULL);
	*/
}

//TODO: the working directory and the executable need to be read from a configuration file
static void runSimulation(uiButton *b, void *data) {


    //TODO: config auto save on run. Make a windows version of this
	#ifdef linux
    if(uiSourceViewGetModified(configText)) {
        //uiMsgBoxError(w, "Invalid file!", "Error parsing ini file!");
        int response = uiMsgBoxConfirmCancel(w, "File not saved!", "The config file was modified. Do you want to save before running?");

        if(response == uiReturnValueCancel) {            
            return;
        }
        else {            
            saveConfigFile();
        }

    }
	#endif

    start_monodomain_exec();
    uiControlDisable(uiControl(btnRun));
}

void onConfigChanged(uiSourceView *e, void *data) {

//TODO: make a windows version of this function

#ifdef linux
    int modified = uiSourceViewGetModified(e);
    int can_undo = uiSourceViewCanUndo(e);

    printf("MOD: %d, CAN: %d\n", modified, can_undo);

    if(modified) {
        uiControlEnable(uiControl(btnSave));
    }
    else {
        uiControlDisable(uiControl(btnSave));
    }
#endif

}

int main()  {

    uiInitOptions o;

    memset(&o, 0, sizeof (uiInitOptions));
    if (uiInit(&o) != NULL)
        abort();

    options = new_user_options ();

    w = uiNewWindow("MonoAlg3D GUI", 1024, 768, 0);
    uiWindowSetMargined(w, 1);

    horizontalButtonBox = uiNewHorizontalBox();
    uiBoxSetPadded(horizontalButtonBox, 1);
    uiWindowSetChild(w, uiControl(horizontalButtonBox));

    horizontalTextBox = uiNewHorizontalBox();
    uiBoxSetPadded(horizontalTextBox, 1);
    uiWindowSetChild(w, uiControl(horizontalTextBox));

    verticalBox = uiNewVerticalBox();
    uiBoxSetPadded(verticalBox, 1);
    uiWindowSetChild(w, uiControl(verticalBox));

	#ifdef linux
    configText = uiNewSourceView();
    uiSourceViewSetReadOnly(configText, true);
    uiSourceViewOnChanged(configText, onConfigChanged, NULL);
	#endif
	
	#ifdef _WIN32
	configText = uiNewMultilineEntry();
    uiMultilineEntrySetReadOnly(configText, 1);
	#endif


    runText = uiNewMultilineEntry();
    uiMultilineEntrySetReadOnly(runText, 1);

    btnSave = uiNewButton("Save");
    uiControlDisable(uiControl(btnSave));
    uiButtonOnClicked(btnSave, saveConfigFileCB, NULL);

    btnConfiguration = uiNewButton("Open configuration");
    uiButtonOnClicked(btnConfiguration, onOpenFileClicked, NULL);

    btnRun = uiNewButton("Run simulation");
    uiButtonOnClicked(btnRun, runSimulation, NULL);
    uiControlDisable(uiControl(btnRun));

    pbar = uiNewProgressBar();

    uiBoxAppend(horizontalTextBox, uiControl(configText), 1);
    uiBoxAppend(horizontalTextBox, uiControl(runText), 1);

    uiBoxAppend(horizontalButtonBox, uiControl(btnConfiguration), 0);
    uiBoxAppend(horizontalButtonBox, uiControl(btnSave), 0);
    uiBoxAppend(horizontalButtonBox, uiControl(btnRun), 0);
    uiBoxAppend(horizontalButtonBox, uiControl(pbar), 1);

    uiBoxAppend(verticalBox, uiControl(horizontalTextBox), 1);
    uiBoxAppend(verticalBox, uiControl(horizontalButtonBox), 0);

    uiWindowOnClosing(w, onClosing, NULL);
    uiControlShow(uiControl(w));
    uiMain();

    return 0;
}
