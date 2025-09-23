//
// Created by sachetto on 18/10/17.
//

#include "file_utils.h"
#include <time.h>
#define STB_DS_IMPLEMENTATION
#include "../3dparty/stb_ds.h"

#include "../3dparty/sds/sds.h"
#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

float *read_timesteps_from_case_file(sds case_file_path) {
    size_t case_file_size;
    float *timesteps = NULL;

    char *content = read_entire_file_with_mmap(case_file_path, &case_file_size);
    char *tmp = content;

    char *word = NULL;

    while(*tmp) {

        if(*tmp != '\n' && *tmp != ':') {
            arrpush(word, *tmp);
        } else {
            arrpush(word, 0);

            if(strcmp(word, "time values") == 0) {
                while(!isdigit(*tmp))
                    tmp++;
                break;
            } else {
                arrsetlen(word, 0);
            }
        }

        tmp++;
    }

    while(*tmp) {
        double time_val = strtod(tmp, &tmp);
        arrput(timesteps, (float)time_val);
        while(isspace(*tmp))
            tmp++;
    }

    munmap(content, case_file_size);
    return timesteps;
}

char *get_timestamped_dir_name(char *dir_name) {

    char *rc = NULL;

    time_t rawtime = time(0);
    struct tm *now = localtime(&rawtime);

    if(rawtime != -1) {
        char timestamp[32];
        strftime(timestamp, 32, "%d_%m_%y_%H_%M_%S", now);
        size_t dir_name_size = strlen(dir_name);
        rc = (char *)malloc(dir_name_size + strlen(timestamp) + 2);
        if(dir_name[dir_name_size - 1] == '/') {
            sprintf(rc, "%.*s_%s", (int)dir_name_size - 1, dir_name, timestamp);
        } else {
            sprintf(rc, "%s_%s", dir_name, timestamp);
        }
    }

    return rc;
}

int get_step_from_filename(char *filename) {
    char *ia = filename;

    int int_a = 0;

    for(; *ia; ia++) {
        if(isdigit(*ia))
            int_a = int_a * 10 + *ia - '0';
    }

    return int_a;
}

FILE *open_file_or_exit(char *filename, char *mode) {

    FILE *f = fopen(filename, mode);

    if(!f) {
        fprintf(stderr, "Error! Could not open %s!\n", filename);
        exit(EXIT_FAILURE);
    }

    return f;
}

char *get_current_directory() {

    char *buf = NULL;
    size_t size = pathconf(".", _PC_PATH_MAX);

    if((buf = (char *)malloc(size)) != NULL)
        buf = getcwd(buf, size);

    buf = strcat(buf, "/");
    return buf;
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename)
        return NULL;
    return dot + 1;
}

char *get_dir_from_path(const char *path) {
    const char *last_slash = NULL;
    char *parent = NULL;
    last_slash = strrchr(path, '/');

    if(last_slash == NULL)
        parent = strdup(".");
    else
        parent = strndup(path, last_slash - path + 1);

    return parent;
}

char *get_file_from_path(const char *path) {
    const char *last_slash = NULL;
    char *file = NULL;
    last_slash = strrchr(path, '/');

    if(last_slash) {
        file = strndup(last_slash + 1, path - last_slash + 1);
        return file;
    } else {
        return strdup(path);
    }
}

char *get_filename_without_ext(const char *filename) {
    const char *last_dot = NULL;
    char *file = NULL;
    last_dot = strrchr(filename, '.');

    if(!last_dot || last_dot == filename)
        return (char *)filename;

    file = strndup(filename, strlen(filename) - strlen(last_dot));
    return file;
}

int cp_file(const char *to, const char *from) {
    int fd_to, fd_from;
    char buf[4096];
    size_t nread;
    int saved_errno;

    fd_from = open(from, O_RDONLY);
    if(fd_from < 0)
        return -1;

    fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if(fd_to < 0)
        goto out_error;

    while(nread = read(fd_from, buf, sizeof buf), nread > 0) {
        char *out_ptr = buf;
        size_t nwritten;

        do {
            nwritten = write(fd_to, out_ptr, nread);

            if(nwritten >= 0) {
                nread -= nwritten;
                out_ptr += nwritten;
            } else if(errno != EINTR) {
                goto out_error;
            }
        } while(nread > 0);
    }

    if(nread == 0) {
        if(close(fd_to) < 0) {
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
    if(fd_to >= 0)
        close(fd_to);

    errno = saved_errno;
    return -1;
}

char *read_entire_file_with_mmap(const char *filename, size_t *size) {

    char *f;

    if(!filename)
        return NULL;

    struct stat s;
    int fd = open(filename, O_RDONLY);

    if(fd == -1) {
        return NULL;
    }

    fstat(fd, &s);
    if(!S_ISREG(s.st_mode)) {
        close(fd);
        return NULL;
    }

    *size = s.st_size;

    size_t to_page_size = *size;

    int pagesize = getpagesize();
    to_page_size += pagesize - (to_page_size % pagesize);

    f = (char *)mmap(0, to_page_size, PROT_READ, MAP_PRIVATE, fd, 0);

    if(f == NULL)
        return NULL;

    close(fd);

    return f;
}

char *read_entire_file(const char *filename, size_t *size) {

    FILE *infile;
    char *buffer;
    size_t numbytes;

    if(!filename)
        return NULL;

    infile = fopen(filename, "r");

    if(infile == NULL)
        return NULL;

    fseek(infile, 0L, SEEK_END);
    numbytes = ftell(infile);

    fseek(infile, 0L, SEEK_SET);

    buffer = (char *)malloc(numbytes * sizeof(char));

    if(buffer == NULL) {
        fclose(infile);
        return NULL;
    }

    fread(buffer, sizeof(char), numbytes, infile);
    fclose(infile);

    *size = numbytes;

    return buffer;
}

string_array read_lines(const char *filename) {

    string_array lines = NULL;

    size_t len = 0;

    FILE *fp;

    fp = fopen(filename, "r");

    if(fp == NULL) {
        fprintf(stderr, "Error reading file %s\n", filename);
        return NULL;
    }

    char *line = NULL;
    while((getline(&line, &len, fp)) != -1) {
        line[strlen(line) - 1] = '\0';
        arrput(lines, strdup(line));
    }

    free(line);
    fclose(fp);

    return lines;
}

/* qsort C-sds comparison function */
static int cstring_cmp(const void *a, const void *b) {
    char *ia = *((char **)a);
    char *ib = *((char **)b);

    int int_a = 0;
    int int_b = 0;

    for(; *ia; ia++) {
        if(isdigit(*ia))
            int_a = int_a * 10 + *ia - '0';
    }

    for(; *ib; ib++) {
        if(isdigit(*ib))
            int_b = int_b * 10 + *ib - '0';
    }

    return int_a - int_b;
    /* strcmp functions works exactly as expected from
    comparison function */
}

// TODO: maybe return the full path of the file?
// We only return the file names here, not the full path!!
string_array list_files_from_dir(const char *dir, const char *prefix, const char *extension, string_array ignore_extensions, bool sort) {

    DIR *dp;

    string_array files = NULL;

    struct dirent *dirp;

    if((dp = opendir(dir)) == NULL) {
        fprintf(stderr, "Error opening %s\n", dir);
        return NULL;
    }

    while((dirp = readdir(dp)) != NULL) {

        if(dirp->d_type != DT_REG) {
            continue;
        }

        char *file_name = strdup(dirp->d_name);

        if(ignore_extensions) {

            int n = arrlen(ignore_extensions);
            bool ignore = false;

            for(int i = 0; i < n; i++) {
                if(FILE_EXTENSION_EQUALS(get_filename_ext(file_name), ignore_extensions[i])) {
                    ignore = true;
                    break;
                }
            }

            if(ignore)
                continue;
        }

        if(extension) {
            const char *file_ext = get_filename_ext(file_name);

            if(!FILE_EXTENSION_EQUALS(file_ext, extension)) {
                continue;
            }
        }

        if(prefix) {
            if(strncmp(prefix, file_name, strlen(prefix)) != 0) {
                continue;
            }
        }

        arrput(files, file_name);
    }

    if(files) {
        if(sort) {
            qsort(files, arrlen(files), sizeof(char *), cstring_cmp);
        }
    }

    closedir(dp);
    return files;
}

bool file_exists(const char *path) {

    if(access(path, F_OK) != -1) {
        // file exists
        return true;
    } else {
        // file doesn't exist
        return false;
    }
}

bool dir_exists(const char *path) {
    struct stat info;

    if(stat(path, &info) != 0) {
        return false;
    } else if(info.st_mode & S_IFDIR) {
        return true;
    }

    return false;
}

void free_path_information(struct path_information *input_info) {
    free(input_info->dir_name);
    free(input_info->filename_without_extension);
    free(input_info->file_extension);
}

void get_path_information(const char *path, struct path_information *input_info) {
    struct stat info;

    input_info->dir_name = NULL;
    input_info->filename_without_extension = NULL;
    input_info->file_extension = NULL;

    input_info->exists = false;

    if(stat(path, &info) != 0) {
        return;
    }

    input_info->exists = true;
    input_info->is_dir = info.st_mode & S_IFDIR;
    input_info->is_file = !(input_info->is_dir);

    if(input_info->is_file) {
        input_info->file_extension = strdup(get_filename_ext(path));
        char *file_name = get_file_from_path(path);
        input_info->filename_without_extension = get_filename_without_ext(file_name);
        input_info->dir_name = get_dir_from_path(path);
        free(file_name);
    } else {
        input_info->dir_name = strdup(path);
    }
}

int remove_directory(const char *path) {
    DIR *d = opendir(path);
    size_t path_len = strlen(path);
    int r = -1;

    if(d) {
        struct dirent *p;

        r = 0;

        while(!r && (p = readdir(d))) {
            int r2 = -1;
            char *buf;
            size_t len;

            /* Skip the names "." and ".." as we don't want to recurse on them. */
            if(!strcmp(p->d_name, ".") || !strcmp(p->d_name, "..")) {
                continue;
            }

            len = path_len + strlen(p->d_name) + 2;
            buf = (char *)malloc(len);

            if(buf) {
                struct stat statbuf;

                snprintf(buf, len, "%s/%s", path, p->d_name);

                if(!stat(buf, &statbuf)) {
                    if(S_ISDIR(statbuf.st_mode)) {
                        r2 = remove_directory(buf);
                    } else {
                        r2 = unlink(buf);
                    }
                }

                free(buf);
            }

            r = r2;
        }

        closedir(d);
    }

    if(!r) {
        r = rmdir(path);
    }

    return r;
}

void create_dir(char *out_dir) {

    if(dir_exists(out_dir))
        return;

    size_t out_dir_len = strlen(out_dir);

    char *new_dir = (char *)malloc(out_dir_len + 2);

    memcpy(new_dir, out_dir, out_dir_len + 1);

    if(new_dir[out_dir_len] != '/') {
        new_dir[out_dir_len] = '/';
        new_dir[out_dir_len + 1] = '\0';
    }

    int start = 0;

    if(new_dir[0] == '/') {
        start++;
    }

    char *slash = strchr(new_dir + start, '/');

    while(slash) {
        size_t dirname_size = slash - new_dir;
        char *dir_to_create = (char *)malloc(dirname_size + 1);
        memcpy(dir_to_create, new_dir, dirname_size);
        dir_to_create[dirname_size] = '\0';

        if(!dir_exists(dir_to_create)) {

            printf("%s does not exist! Creating!\n", dir_to_create);

            if(mkdir(dir_to_create, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) {
                if(!dir_exists(dir_to_create)) { // HACK: this can avoid MPI or batch simulation problems...
                    fprintf(stderr, "Error creating directory %s Exiting!\n", dir_to_create);
                    exit(EXIT_FAILURE);
                }
            }
        }

        slash = strchr(slash + 1, '/');
        free(dir_to_create);
    }

    free(new_dir);
}

/*
 * Base64 encoding/decoding (RFC1341)
 * Copyright (c) 2005-2011, Jouni Malinen <j@w1.fi>
 *
 * This software may be distributed under the terms of the BSD license.
 */
static const unsigned char base64_table[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/**
 * base64_encode - Base64 encode
 * @src: Data to be encoded
 * @len: Length of the data to be encoded
 * @out_len: Pointer to output length variable, or %NULL if not used
 * Returns: Allocated buffer of out_len bytes of encoded data,
 * or %NULL on failure
 *
 * Caller is responsible for freeing the returned buffer. Returned buffer is
 * nul terminated to make it easier to use as a C sds. The nul terminator is
 * not included in out_len.
 */
unsigned char *base64_encode(const unsigned char *src, size_t len, size_t *out_len) {
    unsigned char *out, *pos;
    const unsigned char *end, *in;
    size_t olen;
    int line_len;

    olen = len * 4 / 3 + 4; /* 3-byte blocks to 4-byte */
    olen += olen / 72;      /* line feeds */
    olen++;                 /* nul termination */
    if(olen < len)
        return NULL; /* integer overflow */
    out = (unsigned char *)malloc(olen);
    if(out == NULL)
        return NULL;

    end = src + len;
    in = src;
    pos = out;
    line_len = 0;
    while(end - in >= 3) {
        *pos++ = base64_table[in[0] >> 2];
        *pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
        *pos++ = base64_table[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
        *pos++ = base64_table[in[2] & 0x3f];
        in += 3;
        line_len += 4;
        if(line_len >= 72) {
            *pos++ = '\n';
            line_len = 0;
        }
    }

    if(end - in) {
        *pos++ = base64_table[in[0] >> 2];
        if(end - in == 1) {
            *pos++ = base64_table[(in[0] & 0x03) << 4];
            *pos++ = '=';
        } else {
            *pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
            *pos++ = base64_table[(in[1] & 0x0f) << 2];
        }
        *pos++ = '=';
        line_len += 4;
    }

    if(line_len)
        *pos++ = '\n';

    *pos = '\0';
    if(out_len)
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
size_t base64_decode(unsigned char *out, const char *src, size_t len, size_t *bytes_read) {
    unsigned char dtable[256], *pos, block[4];
    size_t i, count;
    int pad = 0;

    memset(dtable, 0x80, 256);
    for(i = 0; i < sizeof(base64_table) - 1; i++)
        dtable[base64_table[i]] = (unsigned char)i;
    dtable['='] = 0;

    count = 0;
    for(i = 0; i < len; i++) {
        if(dtable[src[i]] != 0x80)
            count++;
    }

    if(count == 0 || count % 4)
        return 0;

    pos = out;
    if(out == NULL)
        return 0;

    count = 0;
    for(i = 0; i < len; i++) {
        unsigned char tmp = dtable[src[i]];
        if(tmp == 0x80)
            continue;

        if(src[i] == '=')
            pad++;
        block[count] = tmp;
        count++;
        if(count == 4) {
            *pos++ = (block[0] << 2) | (block[1] >> 4);
            *pos++ = (block[1] << 4) | (block[2] >> 2);
            *pos++ = (block[2] << 6) | block[3];
            count = 0;
            if(pad) {
                if(pad == 1)
                    pos--;
                else if(pad == 2)
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

bool check_simulation_completed(char *simulation_dir) {

    if(!dir_exists(simulation_dir))
        return false;

    sds output_file_name = sdscatprintf(sdsempty(), "%s/outputlog.txt", simulation_dir);

    char *word = NULL;

    size_t file_size = 0;
    char *outputlog_content = read_entire_file(output_file_name, &file_size);

    if(outputlog_content == NULL) {
        sdsfree(output_file_name);
        return false;
    }

    for(size_t c = 0; c < file_size; c++) {
        char l = outputlog_content[c];
        if(!isspace(l)) {
            arrput(word, l);
        } else {
            arrput(word, '\0');
            if(strcmp(word, "Resolution") == 0) {
                arrfree(word);
                return true;
            }
            arrfree(word);
            word = NULL;
        }
    }

    return false;
}

real_cpu **read_octave_mat_file_to_array(FILE *matrix_file, uint64_t *num_lines, uint64_t *nnz) {
    const char *sep = " ";
    char *line_a = NULL;
    size_t len;
    int count;

    assert(matrix_file);

    do {
        getline(&line_a, &len, matrix_file);
        sds *tmp = sdssplitlen(line_a, (int)strlen(line_a), sep, (int)strlen(sep), &count);
        if(count) {
            if(strcmp(tmp[1], "columns:") == 0) {
                (*num_lines) = strtol(tmp[2], NULL, 10);
            }
            if(strcmp(tmp[1], "nnz:") == 0) {
                (*nnz) = strtol(tmp[2], NULL, 10);
            }
        }
        sdsfreesplitres(tmp, count);
    } while((line_a)[0] == '#');

    real_cpu **matrix = (real_cpu **)malloc(*num_lines * sizeof(real_cpu *));

    for(uint64_t i = 0; i < *num_lines; i++) {
        matrix[i] = (real_cpu *)calloc(*num_lines, sizeof(real_cpu));
    }

    uint64_t item_count = 0;
    uint64_t m_line, m_column;
    real_cpu m_value;

    while(item_count < *nnz) {

        sds *tmp = sdssplitlen(line_a, (int)strlen(line_a), sep, (int)strlen(sep), &count);
        if(tmp[0][0] != '\n') {
            m_line = strtol(tmp[0], NULL, 10);
            m_column = strtol(tmp[1], NULL, 10);
            m_value = (real_cpu)strtod(tmp[2], NULL);

            matrix[m_line - 1][m_column - 1] = m_value;
        }
        sdsfreesplitres(tmp, count);

        item_count++;
        getline(&line_a, &len, matrix_file);
    }

    if(line_a)
        free(line_a);

    return matrix;
}

real_cpu *read_octave_vector_file_to_array(FILE *vec_file, uint64_t *num_lines) {

    ssize_t read;
    size_t len;
    char *line_b = NULL;
    int count;
    char *sep = " ";

    do {
        read = getline(&line_b, &len, vec_file);
        sds *tmp = sdssplitlen(line_b, (int)strlen(line_b), sep, (int)strlen(sep), &count);
        if(count) {
            if(strcmp(tmp[1], "rows:") == 0) {
                (*num_lines) = strtol(tmp[2], NULL, 10);
            }
        }
        sdsfreesplitres(tmp, count);
    } while((line_b)[0] == '#');

    real_cpu *vector = (real_cpu *)malloc(*num_lines * sizeof(real_cpu));

    int item_count = 0;
    while((item_count < *num_lines) && read) {
        sds *tmp = sdssplitlen(line_b, (int)strlen(line_b), sep, (int)strlen(sep), &count);

        if(tmp[0][0] != '\n') {
            vector[item_count] = (real_cpu)strtod(tmp[1], NULL);
        }

        sdsfreesplitres(tmp, count);

        item_count++;
        read = getline(&line_b, &len, vec_file);
    }

    if(line_b)
        free(line_b);

    return vector;
}
