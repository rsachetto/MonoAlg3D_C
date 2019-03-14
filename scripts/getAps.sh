#!/usr/bin/env bash
AP_DIR=$1
AP_PREFIX=$2
AP_LINE=$3
AP_OUT=$4
CELL_ID=$5

if [[ "$#" -ne 5 ]]; then
    echo "-------------------------------------------------------------------"
    echo "Usage:> $0 <AP_DIR> <AP_PREFIX> <AP_LINE> <AP_OUT> <CELL_ID>"
    echo "-------------------------------------------------------------------"
    echo "<AP_DIR> = Directory where the results of the simulation are stored"
    echo "<AP_PREFIX> = Prefix of the results files"
    echo "<AP_LINE> = Line where the action potential information is stored"
    echo "<AP_OUT> = Name for the output text file"
    echo "<CELL_ID> = Index of the cell to be plotted"
    echo "-------------------------------------------------------------------"
    echo "!!! This script only works if the output file is NOT saved in   !!!"
    echo "!!!                     binary format                           !!!"
    echo "-------------------------------------------------------------------"
    echo "<AP_LINE> guide:"
    echo "  VTP/VTU --> Line 6"
    echo "  VTK     --> Last line"
    echo "-------------------------------------------------------------------"
    exit 1
fi

for i in `ls -1v ${AP_DIR}/${AP_PREFIX}*`; do sed -n "${AP_LINE}p" $i | awk -v "col=${CELL_ID}" -F ' ' '{print $col}'  ; done > ${AP_OUT}

