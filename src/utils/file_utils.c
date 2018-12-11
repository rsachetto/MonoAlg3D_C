//
// Created by sachetto on 18/10/17.
//

#include <stdarg.h>

#include "file_utils.h"
#include <stdio.h>
#include <fcntl.h>

#include <errno.h>

#ifdef _WIN32
#include <io.h>
#define read _read
#endif

#ifdef linux
#include <unistd.h>
#include <stdlib.h>

#endif

static FILE *logfile = NULL;

void print_to_stdout_and_file(char const *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    fflush(stdout);
    va_end(ap);
    va_start(ap, fmt);
    if(logfile) {
        vfprintf(logfile, fmt, ap);
        fflush(logfile);
    }
    va_end(ap);
}

void open_logfile(const char *path) {

#ifdef _WIN32
    fopen_s(&logfile, path, "w");
#else
    logfile = fopen(path, "w");
#endif

    if(logfile == NULL) {
        fprintf(stderr, "Error opening %s, printing output only in the sdtout (Terminal)\n", path);
    }
    else {
        printf("Log will be saved in %s\n", path);
    }
}

void close_logfile() {
    if(logfile) fclose(logfile);
}


int cp_file(const char *to, const char *from)
{
    int fd_to, fd_from;
    char buf[4096];
    int nread;
    int saved_errno;

    fd_from = open(from, O_RDONLY);
    if (fd_from < 0)
        return -1;

    fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd_to < 0)
        goto out_error;

    while (nread = read(fd_from, buf, sizeof buf), nread > 0)
    {
        char *out_ptr = buf;
        int  nwritten;

        do {
            nwritten = write(fd_to, out_ptr, nread);

            if (nwritten >= 0)
            {
                nread -= nwritten;
                out_ptr += nwritten;
            }
            else if (errno != EINTR)
            {
                goto out_error;
            }
        } while (nread > 0);
    }

    if (nread == 0)
    {
        if (close(fd_to) < 0)
        {
            fd_to = -1;
            goto out_error;
        }
        close(fd_from);

        /* Success! */
        return 0;
    }

    out_error:
    saved_errno = errno;

    close(fd_from);
    if (fd_to >= 0)
        close(fd_to);

    errno = saved_errno;
    return -1;
}

char *read_entire_file(char *filename, long *size) {

    FILE    *infile;
    char    *buffer;
    long    numbytes;

    if(!filename) return NULL;

/* open an existing file for reading */
    infile = fopen(filename, "r");

/* quit if the file does not exist */
    if(infile == NULL)
        return NULL;

/* Get the number of bytes */
    fseek(infile, 0L, SEEK_END);
    numbytes = ftell(infile);

/* reset the file position indicator to
the beginning of the file */
    fseek(infile, 0L, SEEK_SET);

/* grab sufficient memory for the
buffer to hold the text */
    buffer = (char*)malloc(numbytes*sizeof(char));

/* memory error */
    if(buffer == NULL)
        return NULL;

/* copy all the text into the buffer */
    fread(buffer, sizeof(char), numbytes, infile);
    fclose(infile);

    *size = numbytes;

    return buffer;
}