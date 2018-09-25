// 6 december 2015
#include "uipriv_unix.h"

struct uiSourceView {
	uiUnixControl c;
	GtkWidget *widget;
	GtkContainer *scontainer;
	GtkScrolledWindow *sw;
	GtkWidget *textviewWidget;
	GtkSourceView *sourceView;
	GtkSourceBuffer *sourceBuffer;
	GtkSourceLanguageManager *lm;
	void (*onChanged)(uiSourceView *, void *);
	void *onChangedData;
	gulong onChangedSignal;
};

uiUnixControlAllDefaults(uiSourceView)

static void onChanged(GtkSourceBuffer *textbuf, gpointer data)
{
	uiSourceView *e = uiSourceView(data);
    gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(e->sourceBuffer), TRUE);
	(*(e->onChanged))(e, e->onChangedData);
}

static void defaultOnChanged(uiSourceView *e, void *data)
{
	// do nothing
}

static void saveCallback(GtkSourceFileSaver *saver,
                         GAsyncResult       *result,
                         GtkSourceBuffer    *textbuf) {


    GError *error = NULL;
    gtk_source_file_saver_save_finish (saver, result, &error);

    if (error != NULL)
    {
        g_debug ("File saving error: %s", error->message);
    }

    else {
        gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(textbuf), FALSE);
    }


    if (error != NULL)
    {
        g_error_free (error);
    }


}

char *uiSourceViewText(uiSourceView *e)
{
	GtkTextIter start, end;
	char *tret, *out;
	gtk_text_buffer_get_start_iter(GTK_TEXT_BUFFER(e->sourceBuffer), &start);
	gtk_text_buffer_get_end_iter(GTK_TEXT_BUFFER(e->sourceBuffer), &end);
	tret = gtk_text_buffer_get_text(GTK_TEXT_BUFFER(e->sourceBuffer), &start, &end, TRUE);
	// theoretically we could just return tret because uiUnixStrdupText() is just g_strdup(), but if that ever changes we can't, so let's do it this way to be safe
	out = uiUnixStrdupText(tret);
	g_free(tret);
	return out;
}

int uiSourceViewGetModified(uiSourceView *e) {
    return gtk_text_buffer_get_modified(GTK_TEXT_BUFFER(e->sourceBuffer));
}

void uiSourceViewSetText(uiSourceView *e, const char *text, char *mime_type)
{

    GtkSourceLanguage *language = NULL;
	language = gtk_source_language_manager_guess_language(e->lm, NULL, g_content_type_from_mime_type(mime_type));
	gtk_source_buffer_set_language (e->sourceBuffer, language);

	// we need to inhibit sending of ::changed because this WILL send a ::changed otherwise
	g_signal_handler_block(e->sourceBuffer, e->onChangedSignal);

	gtk_source_buffer_begin_not_undoable_action (e->sourceBuffer);
	gtk_text_buffer_set_text (GTK_TEXT_BUFFER (e->sourceBuffer), text, -1);
    gtk_source_buffer_end_not_undoable_action (e->sourceBuffer);

    g_signal_handler_unblock(e->sourceBuffer, e->onChangedSignal);
}

void uiSourceViewOnChanged(uiSourceView *e, void (*f)(uiSourceView *e, void *data), void *data)
{
	e->onChanged = f;
	e->onChangedData = data;
}

static uiSourceView *finishSourceView(GtkPolicyType hpolicy)
{
	uiSourceView *e;

	uiUnixNewControl(uiSourceView, e);

	e->widget = gtk_scrolled_window_new(NULL, NULL);
	e->scontainer = GTK_CONTAINER(e->widget);
	e->sw = GTK_SCROLLED_WINDOW(e->widget);
	gtk_scrolled_window_set_policy(e->sw,
		hpolicy,
		GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_set_shadow_type(e->sw, GTK_SHADOW_IN);


	e->sourceBuffer = GTK_SOURCE_BUFFER (gtk_source_buffer_new (NULL));

	e->textviewWidget = gtk_source_view_new_with_buffer(e->sourceBuffer);
	e->sourceView = GTK_SOURCE_VIEW(e->textviewWidget);

	e->lm = gtk_source_language_manager_new();
	gtk_container_add(e->scontainer, e->textviewWidget);
	// and make the text view visible; only the scrolled window's visibility is controlled by libui
	gtk_widget_show(e->textviewWidget);

	e->onChangedSignal = g_signal_connect(e->sourceBuffer, "changed", G_CALLBACK(onChanged), e);
	uiSourceViewOnChanged(e, defaultOnChanged, NULL);

	return e;
}

void uiSourceViewSetReadOnly(uiSourceView *e, int readonly)
{
    gboolean editable;

    editable = TRUE;
    if (readonly)
        editable = FALSE;
    gtk_text_view_set_editable(GTK_TEXT_VIEW(e->sourceView), editable);
}

int uiSourceViewCanUndo(uiSourceView *e)
{
    return gtk_source_buffer_can_undo(e->sourceBuffer);
}


//TODO: handle error better
void uiSourceViewSaveSource(uiSourceView *e, const char * filename) {

    GFile *file = g_file_new_for_path(filename);
    GtkSourceFile *sourceFile = gtk_source_file_new();

    gtk_source_file_set_location(sourceFile, file);

    GtkSourceFileSaver *fileSaver = gtk_source_file_saver_new(e->sourceBuffer, sourceFile);

    gtk_source_file_saver_save_async(fileSaver, G_PRIORITY_DEFAULT, NULL, NULL, NULL, NULL, (GAsyncReadyCallback) saveCallback, e->sourceBuffer);

}

uiSourceView *uiNewSourceView(void)
{
	return finishSourceView(GTK_POLICY_NEVER);
}

