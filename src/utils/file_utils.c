//
// Created by sachetto on 18/10/17.
//

#include <stdarg.h>

#include "file_utils.h"
#include "../string/sds.h"
#include <stdio.h>
#include <fcntl.h>
#include <string.h>

#include <errno.h>

#include "../single_file_libraries/stb_ds.h"


#ifdef _WIN32
#include <io.h>
#include <Windows.h>
#include <direct.h>
#include <sys/stat.h>
#define read _read
#endif

#ifdef linux

#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>

#endif

static FILE *logfile = NULL;

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

char *read_entire_file(char *filename, long *size) {

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

#ifdef _WIN32
// if typedef doesn't exist (msvc, blah)
typedef intptr_t ssize_t;

ssize_t getline(char **lineptr, size_t *n, FILE *stream) {
    size_t pos;
    int c;

    if (lineptr == NULL || stream == NULL || n == NULL) {
        errno = EINVAL;
        return -1;
    }

    c = fgetc(stream);
    if (c == EOF) {
        return -1;
    }

    if (*lineptr == NULL) {
        *lineptr = malloc(128);
        if (*lineptr == NULL) {
            return -1;
        }
        *n = 128;
    }

    pos = 0;
    while(c != EOF) {
        if (pos + 1 >= *n) {
            size_t new_size = *n + (*n >> 2);
            if (new_size < 128) {
                new_size = 128;
            }
            char *new_ptr = realloc(*lineptr, new_size);
            if (new_ptr == NULL) {
                return -1;
            }
            *n = new_size;
            *lineptr = new_ptr;
        }

        ((unsigned char *)(*lineptr))[pos ++] = c;
        if (c == '\n') {
            break;
        }
        c = fgetc(stream);
    }

    (*lineptr)[pos] = '\0';
    return pos;
}
#endif

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


#ifndef _WIN32
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
#endif

#ifdef _WIN32
bool dir_exists(const char *path)
{
    DWORD ftyp = GetFileAttributesA(path);
    if (ftyp == INVALID_FILE_ATTRIBUTES)
        return false;  //something is wrong with your path!

    if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
        return true;   // this is a directory!

    return false;    // this is not a directory!
}

int is_dots(const char* str) {
    if(strcmp(str,".") && strcmp(str,"..")) return FALSE;
    return TRUE;
}

int remove_directory(const char *path) {
    HANDLE hFind;    // file handle
    WIN32_FIND_DATA FindFileData;

    char DirPath[MAX_PATH];
    char FileName[MAX_PATH];

    strcpy(DirPath,path);
    strcat(DirPath,"\\*");    // searching all files
    strcpy(FileName,path);
    strcat(FileName,"\\");

    // find the first file
    hFind = FindFirstFile(DirPath,&FindFileData);
    if(hFind == INVALID_HANDLE_VALUE) return FALSE;
    strcpy(DirPath,FileName);

    bool bSearch = true;
    while(bSearch) {    // until we find an entry
        if(FindNextFile(hFind,&FindFileData)) {
            if(is_dots(FindFileData.cFileName)) continue;
            strcat(FileName,FindFileData.cFileName);
            if((FindFileData.dwFileAttributes &
                FILE_ATTRIBUTE_DIRECTORY)) {

                // we have found a directory, recurse
                if(!remove_directory(FileName)) {
                    FindClose(hFind);
                    return FALSE;    // directory couldn't be deleted
                }
                // remove the empty directory
                RemoveDirectory(FileName);
                strcpy(FileName,DirPath);
            }
            else {
                if(FindFileData.dwFileAttributes &
                   FILE_ATTRIBUTE_READONLY)
                    // change read-only file mode
                    _chmod(FileName, _S_IWRITE);
                if(!DeleteFile(FileName)) {    // delete the file
                    FindClose(hFind);
                    return FALSE;
                }
                strcpy(FileName,DirPath);
            }
        }
        else {
            // no more files there
            if(GetLastError() == ERROR_NO_MORE_FILES)
                bSearch = false;
            else {
                // some error occurred; close the handle and return FALSE
                FindClose(hFind);
                return FALSE;
            }

        }

    }
    FindClose(hFind);                  // close the file handle

    return RemoveDirectory(path);     // remove the empty directory

}

#else
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

#endif
//
//void create_dir(const char *out_dir) {
//
//    //TODO: check for windows dir separators
//    int dirs_count;
//
//    sds *all_dirs = sdssplit(out_dir, "/", &dirs_count);
//    sds new_dir = sdsempty();
//
//    for(int d = 0; d < dirs_count; d++) {
//
//        new_dir = sdscat(new_dir, all_dirs[d]);
//        new_dir = sdscat(new_dir, "/");
//
//        if (!dir_exists (new_dir)) {
//
//            printf ("%s does not exist! Creating!\n", new_dir);
//#if defined _MSC_VER
//            if (_mkdir(out_dir) == -1)
//#else
//            if (mkdir(new_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
//#endif
//            {
//
//                fprintf (stderr, "Error creating directory %s Exiting!\n", new_dir);
//                exit (EXIT_FAILURE);
//            }
//        }
//    }
//
//    sdsfreesplitres(all_dirs, dirs_count);
//    sdsfree(new_dir);
//
//}

void fixpath(char *path)
{
    for(; *path; ++path)
        if (*path == '\\')
            *path = '/';
}

void create_dir(char *out_dir) {

    //TODO: check why this is not working in windows. It seems that out_dir is not null terminated.
    //fixpath(out_dir);

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
#if defined _MSC_VER
            if (_mkdir(dir_to_create) == -1)
#else
                if (mkdir(dir_to_create, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
#endif
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


