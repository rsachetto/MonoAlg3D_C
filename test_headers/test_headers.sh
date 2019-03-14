#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    [ -s test_h ] || rm test_h

    error_file=`basename ${line}`

    echo "#include \"$line\" " > test_h.c
    echo "#include <stdio.h> " >> test_h.c
    echo "" >> test_h.c

    echo "int main() {printf(\"OHOHOHO\");}" >> test_h.c

    echo "" >> test_h.c
    cc  -I/opt/cuda/include/ -Werror  test_h.c   -o test_h 2> ${error_file}.txt

    if [ $? -eq 0 ]; then
        echo Header $line compiled sucessfuly!
    else
        echo Failed to compile header $line
    fi


    if [ -s ${error_file}.txt ]; then
        cat test_h.c >>  ${error_file}.txt
    else
         rm ${error_file}.txt
    fi


done < "$1"
