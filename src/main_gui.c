#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <stdbool.h>
#include "ini_parser/ini.h"
#include "config/config_parser.h"
#include "fork_pipe/fork_pipe.h"

#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#endif

#include "libui/ui.h"


#ifdef linux
#include <unistd.h>
#endif

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


static void appendToOutput(void *data) {
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

static void appendToQueue(void *data) {
    uiQueueMain(appendToOutput, data);
}


static int onClosing(uiWindow *w, void *data) {
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

struct ThreadData {
    char *program;
    void (*fn_pointer)(void*);
};


#ifdef _WIN32
DWORD WINAPI start_program_with_tread(LPVOID thread_param) {
#else
void start_program_with_tread(void * thread_param) {
#endif

    struct ThreadData *td = (struct ThreadData*) thread_param;
    run_child_process_and_process_output(td->program, td->fn_pointer);
}


void start_monodomain_exec() {

    //LEAK: maybe free it on the run_child_process_and_process_output function
    struct ThreadData *td = (struct ThreadData *)malloc(sizeof(struct ThreadData));
    td->fn_pointer = appendToQueue;

#ifdef _WIN32
    HANDLE thread;
    DWORD  threadId;

    //TODO: we need to pass the Monoalg path and parameters here
    char program[] = "bin\\MonoAlg3D.exe -c example_configs\\benchmark_config_example.ini";
    td->program = _strdup(program);

    thread = CreateThread(
            NULL,                   // default security attributes
            0,                      // use default stack size
            start_program_with_tread,       // thread function name
            (LPVOID) td,          // argument to thread function
            0,                      // use default creation flags
            &threadId);   // returns the thread identifier

#else
    //TODO: we need to pass the Monoalg path and parameters here
    char **program = (char**) malloc(4*sizeof(char*));

    program[0] = strdup("bin/MonoAlg3D");
    program[1] = strdup("-c");
    program[2] = strdup(config_file_name);
    printf ("%s\n", config_file_name);

    program[3] = NULL;



#endif


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
