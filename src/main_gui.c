#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <ui.h>
#include <glib.h>
#include <assert.h>
#include <stdbool.h>
#include "ini_parser/ini.h"
#include "config/config_parser.h"

uiMultilineEntry *configText, *runText;
uiWindow *w;
uiBox *verticalBox, *horizontalButtonBox, *horizontalTextBox;
uiButton *btnConfiguration, *btnRun;
struct user_options *options;

char *config_file_name = NULL;
char *output_last_sim = NULL;

GIOChannel *channel_stderr, *channel_stdout;

int onClosing(uiWindow *w, void *data)
{
    uiQuit();
    return 1;
}

static void onOpenFileClicked(uiButton *b, void *data)
{

    if(config_file_name) uiFreeText(config_file_name);

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
        uiMultilineEntrySetText(configText, string);
        uiControlEnable(uiControl(btnRun));
    }
    else {
        uiMsgBoxError(w, "Invalid file!", "Error parsing ini file!");
    }
}

static gboolean getOutput (GIOChannel *channel, GIOCondition cond, gpointer data) {
    gchar *str_return;
    gsize length;
    gsize terminator_pos;
    GError *error = NULL;

    if (g_io_channel_read_line (channel, &str_return, &length, &terminator_pos, &error) == G_IO_STATUS_ERROR) {
        g_warning ("Something went wrong");
    }
    if (error != NULL) {
        g_warning (error->message);
    }

    uiMultilineEntryAppend(runText, str_return);

    g_free( str_return );

    return TRUE;
}

static void child_watch_cb (GPid pid, gint status, gpointer user_data) {
    g_message ("Child %" G_PID_FORMAT " exited %s", pid,
               g_spawn_check_exit_status (status, NULL) ? "normally" : "abnormally");

    g_io_channel_shutdown(channel_stderr,TRUE,NULL);
    g_io_channel_shutdown(channel_stdout,TRUE,NULL);
    g_spawn_close_pid (pid);
}

static void runSimulation(uiButton *b, void *data) {

    //TODO: the working directory and the executable need to be read from a configuration file

    char **program = (char**) malloc(4*sizeof(char*));

    program[0] = strdup("bin/MonoAlg3D");
    program[1] = strdup("-c");
    program[2] = strdup(config_file_name);
    printf ("%s\n", config_file_name);

    program[3] = NULL;

    gint child_stdout, child_stderr;
    GPid child_pid;
    g_autoptr(GError) error = NULL;

    // Spawn child process.
    g_spawn_async_with_pipes ("/home/sachetto/Projects/MonoAlg3D_C/", program, NULL, G_SPAWN_DO_NOT_REAP_CHILD, NULL,
                              NULL, &child_pid, NULL, &child_stdout,
                              &child_stderr, &error);
    if (error != NULL)
    {
        g_error ("Spawning child failed: %s", error->message);
    }

    g_child_watch_add (child_pid, child_watch_cb, NULL);


    channel_stderr = g_io_channel_unix_new(child_stderr);
    g_io_add_watch (channel_stderr, G_IO_IN, getOutput, NULL);

    channel_stdout = g_io_channel_unix_new(child_stdout);
    g_io_add_watch (channel_stdout, G_IO_IN, getOutput, NULL);
    //g_io_channel_shutdown(channel,TRUE,NULL);

}

int main(void)  {

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

    configText = uiNewMultilineEntry();
    //uiMultilineEntrySetReadOnly(configText, 1);

    runText = uiNewMultilineEntry();
    uiMultilineEntrySetReadOnly(runText, 1);

    btnConfiguration = uiNewButton("Open configuration");
    uiButtonOnClicked(btnConfiguration, onOpenFileClicked, NULL);

    btnRun = uiNewButton("Run simulation");
    uiButtonOnClicked(btnRun, runSimulation, NULL);
    uiControlDisable(uiControl(btnRun));

    uiBoxAppend(horizontalTextBox, uiControl(configText), 1);
    uiBoxAppend(horizontalTextBox, uiControl(runText), 1);

    uiBoxAppend(horizontalButtonBox, uiControl(btnConfiguration), 0);
    uiBoxAppend(horizontalButtonBox, uiControl(btnRun), 0);

    uiBoxAppend(verticalBox, uiControl(horizontalTextBox), 1);
    uiBoxAppend(verticalBox, uiControl(horizontalButtonBox), 0);

    uiWindowOnClosing(w, onClosing, NULL);
    uiControlShow(uiControl(w));
    uiMain();
    return 0;
}
