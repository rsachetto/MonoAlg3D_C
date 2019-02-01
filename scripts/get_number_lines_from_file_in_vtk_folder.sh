#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "-------------------------------------------------------------------"
    echo "Usage:> $0 <AP_DIR> <AP_PREFIX>"
    echo "-------------------------------------------------------------------"
    echo "<AP_DIR> = Directory where the results of the simulation are stored"
    echo "<AP_PREFIX> = Prefix of the results files"
    echo "-------------------------------------------------------------------"
    echo "!!! This script only works if the output file is NOT saved in   !!!"
    echo "!!!                     binary format                           !!!"
    echo "-------------------------------------------------------------------"
    exit 1
fi

AP_DIR=$1
AP_PREFIX=$2

FIRST_FILE=$(ls -1v ${AP_DIR}/${AP_PREFIX}_it_0_*)
NUMBER_OF_LINES=$(wc $FIRST_FILE | awk -F ' ' '{print $1}')
NUMBER_OF_LINES=$((NUMBER_OF_LINES+1))
echo $NUMBER_OF_LINES

