// 26 june 2015
#include "uipriv_unix.h"

// LONGTERM figure out why, and describe, that this is the desired behavior
// LONGTERM also point out that font and color buttons also work like this

#define windowWindow(w) (GTK_WINDOW(uiControlHandle(uiControl(w))))

static const struct {
	gint gtkreturn;
	uiDialogReturnValues extreturn;
} extResponses[] = {
        { GTK_RESPONSE_OK, uiReturnValueOk},
		{ GTK_RESPONSE_CANCEL, uiReturnValueCancel},
		{ GTK_RESPONSE_ACCEPT, uiReturnValueAccept},
        { 0xFFFF, 0 },

};


static gint uiResponseFromExtResponse(gint gtkreturn)
{
    int i;
    for (i = 0; extResponses[i].gtkreturn != 0xFFFF; i++) {
        if (extResponses[i].gtkreturn == gtkreturn) {
            return extResponses[i].extreturn;
        }
    }

    return 0;
}

static char *filedialog(GtkWindow *parent, GtkFileChooserAction mode, const gchar *confirm)
{
	GtkWidget *fcd;
	GtkFileChooser *fc;
	gint response;
	char *filename;

	fcd = gtk_file_chooser_dialog_new(NULL, parent, mode,
		"_Cancel", GTK_RESPONSE_CANCEL,
		confirm, GTK_RESPONSE_ACCEPT,
		NULL);
	fc = GTK_FILE_CHOOSER(fcd);
	gtk_file_chooser_set_local_only(fc, FALSE);
	gtk_file_chooser_set_select_multiple(fc, FALSE);
	gtk_file_chooser_set_show_hidden(fc, TRUE);
	gtk_file_chooser_set_do_overwrite_confirmation(fc, TRUE);
	gtk_file_chooser_set_create_folders(fc, TRUE);
	response = gtk_dialog_run(GTK_DIALOG(fcd));
	if (response != GTK_RESPONSE_ACCEPT) {
		gtk_widget_destroy(fcd);
		return NULL;
	}
	filename = uiUnixStrdupText(gtk_file_chooser_get_filename(fc));
	gtk_widget_destroy(fcd);
	return filename;
}

char *uiOpenFile(uiWindow *parent)
{
	return filedialog(windowWindow(parent), GTK_FILE_CHOOSER_ACTION_OPEN, "_Open");
}

char *uiSaveFile(uiWindow *parent)
{
	return filedialog(windowWindow(parent), GTK_FILE_CHOOSER_ACTION_SAVE, "_Save");
}

static int msgbox(GtkWindow *parent, const char *title, const char *description, GtkMessageType type, GtkButtonsType buttons)
{
	GtkWidget *md;

	md = gtk_message_dialog_new(parent, GTK_DIALOG_MODAL,
		type, buttons,
		"%s", title);
	gtk_message_dialog_format_secondary_text(GTK_MESSAGE_DIALOG(md), "%s", description);
	gint response = gtk_dialog_run(GTK_DIALOG(md));
	gtk_widget_destroy(md);

	return response;
}

void uiMsgBox(uiWindow *parent, const char *title, const char *description)
{
	msgbox(windowWindow(parent), title, description, GTK_MESSAGE_OTHER, GTK_BUTTONS_OK);
}

void uiMsgBoxError(uiWindow *parent, const char *title, const char *description)
{
	msgbox(windowWindow(parent), title, description, GTK_MESSAGE_ERROR, GTK_BUTTONS_OK);
}

int uiMsgBoxConfirmCancel(uiWindow *parent, const char *title, const char *description)
{
	int response = msgbox(windowWindow(parent), title, description, GTK_MESSAGE_QUESTION, GTK_BUTTONS_OK_CANCEL);
	return uiResponseFromExtResponse(response);
}
