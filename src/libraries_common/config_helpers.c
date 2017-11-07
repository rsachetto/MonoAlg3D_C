//
// Created by sachetto on 22/10/17.
//

#include "config_helpers.h"

char *get_char_parameter (struct string_hash *config, const char *parameter) {
    char *config_char = string_hash_search (config, parameter);
    return config_char;
}