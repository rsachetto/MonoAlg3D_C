//
// Created by sachetto on 18/10/17.
//

#include "file_utils.h"

#define STB_DS_IMPLEMENTATION
#include "../single_file_libraries/stb_ds.h"


#include "../string/sds.h"
#include <stdarg.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>


static FILE *logfile = NULL;

char * get_dir_from_path(const char * path) {
    char *last_slash = NULL;
    char *parent = NULL;
    last_slash = strrchr(path, '/');
    parent = strndup(path, last_slash - path + 1);
    return parent;
}

void print_to_stdout_and_file(char const *fmt, ...) {
    va_list ap;

    if (!no_stdout) {
        va_start(ap, fmt);
        vprintf(fmt, ap);
        fflush(stdout);
        va_end(ap);
    }

    va_start(ap, fmt);
    if (logfile) {
        vfprintf(logfile, fmt, ap);
        fflush(logfile);
    }
    va_end(ap);
}

void print_to_stderr_and_file_and_exit(char const *fmt, ...) {
    va_list ap;
            va_start(ap, fmt);
    vprintf(fmt, ap);
    fflush(stderr);
            va_end(ap);
            va_start(ap, fmt);
    if (logfile) {
        vfprintf(logfile, fmt, ap);
        fflush(logfile);
    }
            va_end(ap);
    exit(EXIT_FAILURE);
}

void open_logfile(const char *path) {

#ifdef _WIN32
    fopen_s(&logfile, path, "w");
#else
    logfile = fopen(path, "w");
#endif

    if (logfile == NULL) {
        fprintf(stderr, "Error opening %s, printing output only in the sdtout (Terminal)\n", path);
    } else {
        printf("Log will be saved in %s\n", path);
    }
}

void close_logfile() {
    if (logfile) fclose(logfile);
}


int cp_file(const char *to, const char *from) {
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

    while (nread = read(fd_from, buf, sizeof buf), nread > 0) {
        char *out_ptr = buf;
        int nwritten;

        do {
            nwritten = write(fd_to, out_ptr, nread);

            if (nwritten >= 0) {
                nread -= nwritten;
                out_ptr += nwritten;
            } else if (errno != EINTR) {
                goto out_error;
            }
        } while (nread > 0);
    }

    if (nread == 0) {
        if (close(fd_to) < 0) {
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

char *read_entire_file_with_mmap(const char *filename, size_t *size) {

    char *f;


    if (!filename) return NULL;

    struct stat s;
    int fd = open (filename, O_RDONLY);

    if(fd == -1) {
        return NULL;
    }

    /* Get the size of the file. */
    fstat (fd, & s);
    *size = s.st_size;

    size_t to_page_size = *size;

    int pagesize = getpagesize();
    to_page_size += pagesize - (to_page_size%pagesize);

    f = (char *) mmap (0, to_page_size, PROT_READ, MAP_PRIVATE, fd, 0);

    close(fd);

    if (f == NULL)
        return NULL;

    return f;
}

char *read_entire_file(const char *filename, long *size) {

    FILE *infile;
    char *buffer;
    long numbytes;

    if (!filename) return NULL;

/* open an existing file for reading */
    infile = fopen(filename, "r");

/* quit if the file does not exist */
    if (infile == NULL)
        return NULL;

/* Get the number of bytes */
    fseek(infile, 0L, SEEK_END);
    numbytes = ftell(infile);

/* reset the file position indicator to
the beginning of the file */
    fseek(infile, 0L, SEEK_SET);

/* grab sufficient memory for the
buffer to hold the text */
    buffer = (char *) malloc(numbytes * sizeof(char));

/* memory error */
    if (buffer == NULL)
        return NULL;

/* copy all the text into the buffer */
    fread(buffer, sizeof(char), numbytes, infile);
    fclose(infile);

    *size = numbytes;

    return buffer;
}

string_array read_lines(const char *filename) {

    string_array lines = NULL;


    size_t len = 0;
    ssize_t read;

    FILE *fp;

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Error reading file %s\n", filename);
        return NULL;
    }

    char * line = NULL;
    while ((read = getline(&line, &len, fp)) != -1) {
        line[strlen(line) - 1] = '\0';
        arrput(lines, strdup(line));
    }

    free(line);
    fclose(fp);

    return lines;

}


string_array list_files_from_dir(const char *dir, const char *prefix) {

    DIR *dp;

    string_array files = NULL;

    struct dirent *dirp;

    if ((dp = opendir(dir)) == NULL) {
        fprintf(stderr, "Error opening %s\n", dir);
        return NULL;
    }

    while ((dirp = readdir(dp)) != NULL) {

        char *file_name = strdup(dirp->d_name);

        if (prefix) {

            if (strncmp(prefix, file_name, strlen(prefix)) == 0) {
                arrput(files, file_name);
            }

        } else {
            arrput(files, file_name);
        }
    }

    closedir(dp);
    return files;
}

/* qsort C-string comparison function */
static int cstring_cmp(const void *a, const void *b)
{
    char *ia = *((char **)a);
    char *ib = *((char **)b);

    int int_a = 0;
    int int_b = 0;

    for(; *ia; ia++) {
        if(isdigit(*ia))
            int_a = int_a*10 + *ia - '0';
    }

    for(; *ib; ib++) {
        if(isdigit(*ib))
            int_b = int_b*10 + *ib - '0';
    }

    return int_a - int_b;
    /* strcmp functions works exactly as expected from
    comparison function */
}



string_array list_files_from_dir_sorted(const char *dir, const char *prefix) {

    DIR *dp;

    string_array files = NULL;

    struct dirent *dirp;

    if ((dp = opendir(dir)) == NULL) {
        fprintf(stderr, "Error opening %s\n", dir);
        return NULL;
    }

    while ((dirp = readdir(dp)) != NULL) {

        char *file_name = strdup(dirp->d_name);

        if (prefix) {

            if (strncmp(prefix, file_name, strlen(prefix)) == 0) {
                arrput(files, file_name);
            }

        } else {
            arrput(files, file_name);
        }
    }

    qsort(files, arrlen(files), sizeof(char *), cstring_cmp);

    closedir(dp);
    return files;
}


bool dir_exists(const char *path) {
    struct stat info;

    if(stat( path, &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}

int remove_directory(const char *path)
{
    DIR *d = opendir(path);
    size_t path_len = strlen(path);
    int r = -1;

    if (d)
    {
        struct dirent *p;

        r = 0;

        while (!r && (p=readdir(d)))
        {
            int r2 = -1;
            char *buf;
            size_t len;

            /* Skip the names "." and ".." as we don't want to recurse on them. */
            if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, ".."))
            {
                continue;
            }

            len = path_len + strlen(p->d_name) + 2;
            buf = malloc(len);

            if (buf)
            {
                struct stat statbuf;

                snprintf(buf, len, "%s/%s", path, p->d_name);

                if (!stat(buf, &statbuf))
                {
                    if (S_ISDIR(statbuf.st_mode))
                    {
                        r2 = remove_directory(buf);
                    }
                    else
                    {
                        r2 = unlink(buf);
                    }
                }

                free(buf);
            }

            r = r2;
        }

        closedir(d);
    }

    if (!r)
    {
        r = rmdir(path);
    }

    return r;
}

void fixpath(char *path)
{
    for(; *path; ++path)
        if (*path == '\\')
            *path = '/';
}

void create_dir(char *out_dir) {

    if(dir_exists(out_dir)) return;

    size_t out_dir_len = strlen(out_dir);

    char *new_dir = (char*) malloc(out_dir_len+2);

    memcpy(new_dir, out_dir, out_dir_len+1);

    if(new_dir[out_dir_len] != '/') {
        new_dir[out_dir_len] = '/';
        new_dir[out_dir_len+1] = '\0';
    }

    int start = 0;

    if(new_dir[0] == '/') {
        start++;
    }

    char *slash = strchr(new_dir + start, '/');

    while(slash) {
        size_t dirname_size = slash - new_dir;
        char *dir_to_create = malloc(dirname_size + 1);
        memcpy(dir_to_create, new_dir, dirname_size);
        dir_to_create[dirname_size] = '\0';

        if (!dir_exists (dir_to_create)) {

            printf ("%s does not exist! Creating!\n", dir_to_create);

            if (mkdir(dir_to_create, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
            {
                fprintf (stderr, "Error creating directory %s Exiting!\n", dir_to_create);
                exit (EXIT_FAILURE);
            }
        }

        slash = strchr(slash+1,'/');
        free(dir_to_create);
    }

    free(new_dir);
}


/*
 * Base64 encoding/decoding (RFC1341)
 * Copyright (c) 2005-2011, Jouni Malinen <j@w1.fi>
 *
 * This software may be distributed under the terms of the BSD license.
 * See README for more details.
 */
static const unsigned char base64_table[65] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/**
 * base64_encode - Base64 encode
 * @src: Data to be encoded
 * @len: Length of the data to be encoded
 * @out_len: Pointer to output length variable, or %NULL if not used
 * Returns: Allocated buffer of out_len bytes of encoded data,
 * or %NULL on failure
 *
 * Caller is responsible for freeing the returned buffer. Returned buffer is
 * nul terminated to make it easier to use as a C string. The nul terminator is
 * not included in out_len.
 */
unsigned char * base64_encode(const unsigned char *src, size_t len,
                              size_t *out_len)
{
    unsigned char *out, *pos;
    const unsigned char *end, *in;
    size_t olen;
    int line_len;

    olen = len * 4 / 3 + 4; /* 3-byte blocks to 4-byte */
    olen += olen / 72; /* line feeds */
    olen++; /* nul termination */
    if (olen < len)
        return NULL; /* integer overflow */
    out = malloc(olen);
    if (out == NULL)
        return NULL;

    end = src + len;
    in = src;
    pos = out;
    line_len = 0;
    while (end - in >= 3) {
        *pos++ = base64_table[in[0] >> 2];
        *pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
        *pos++ = base64_table[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
        *pos++ = base64_table[in[2] & 0x3f];
        in += 3;
        line_len += 4;
        if (line_len >= 72) {
            *pos++ = '\n';
            line_len = 0;
        }
    }

    if (end - in) {
        *pos++ = base64_table[in[0] >> 2];
        if (end - in == 1) {
            *pos++ = base64_table[(in[0] & 0x03) << 4];
            *pos++ = '=';
        } else {
            *pos++ = base64_table[((in[0] & 0x03) << 4) |
                                  (in[1] >> 4)];
            *pos++ = base64_table[(in[1] & 0x0f) << 2];
        }
        *pos++ = '=';
        line_len += 4;
    }

    if (line_len)
        *pos++ = '\n';

    *pos = '\0';
    if (out_len)
        *out_len = pos - out;
    return out;
}


/**
 * base64_decode - Base64 decode
 * @src: Data to be decoded
 * @len: Length of the data to be decoded
 * @bytes_read: Pointer to the bytes read on the src stream
 * Returns: out_len bytes of decoded data,
 * or 0 on failure
 *
 * Caller is responsible for allocating the out buffer.
 */
size_t base64_decode(unsigned char *out, const char *src, size_t len, size_t *bytes_read)
{
    unsigned char dtable[256], *pos, block[4], tmp;
    size_t i, count;
    int pad = 0;

    memset(dtable, 0x80, 256);
    for (i = 0; i < sizeof(base64_table) - 1; i++)
        dtable[base64_table[i]] = (unsigned char) i;
    dtable['='] = 0;

    count = 0;
    for (i = 0; i < len; i++) {
        if (dtable[src[i]] != 0x80)
            count++;
    }

    if (count == 0 || count % 4)
        return 0;

    pos = out;
    if (out == NULL)
        return 0;

    count = 0;
    for (i = 0; i < len; i++) {
        tmp = dtable[src[i]];
        if (tmp == 0x80)
            continue;

        if (src[i] == '=')
            pad++;
        block[count] = tmp;
        count++;
        if (count == 4) {
            *pos++ = (block[0] << 2) | (block[1] >> 4);
            *pos++ = (block[1] << 4) | (block[2] >> 2);
            *pos++ = (block[2] << 6) | block[3];
            count = 0;
            if (pad) {
                if (pad == 1)
                    pos--;
                else if (pad == 2)
                    pos -= 2;
                else {
                    return 0;
                }
                break;
            }
        }
    }

    *bytes_read = i + 1;

    return pos - out;
}

