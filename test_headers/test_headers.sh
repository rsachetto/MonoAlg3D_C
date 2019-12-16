#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do

    error_file=$(basename "${line}")

    if cc  -I/opt/cuda/include/ -Werror "${line}" -o tmp.gch > "${error_file}".txt
    then
        echo Header "$line" compiled sucessfuly!
        rm "${error_file}".txt
    else
        echo Failed to compile header "$line"
    fi

    rm tmp.gch

done < "$1"