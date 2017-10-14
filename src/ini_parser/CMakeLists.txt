set(SOURCE_FILES_INI_PARSER
        ini_file_sections.h
        ini.c
        ini.h
        )

ADD_LIBRARY(ini_parser STATIC ${SOURCE_FILES_INI_PARSER})
