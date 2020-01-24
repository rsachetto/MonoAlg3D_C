#!/usr/bin/env bash
APD_DIR=$1
APD_LINE=$2
APD_OUT=$3
CELL_ID=$4

if [[ "$#" -ne 4 ]]; then
    echo "-------------------------------------------------------------------"
    echo "Usage:> $0 <APD_DIR> <APD_LINE> <APD_OUT> <CELL_ID>"
    echo "-------------------------------------------------------------------"
    echo "<APD_DIR> = Directory where the results of the simulation are stored"
    echo "<APD_LINE> = Line where the APD information is stored"
    echo "<APD_OUT> = Name for the output text file"
    echo "<CELL_ID> = Index of the cell to be plotted"
    echo "-------------------------------------------------------------------"
    echo "!!! This script only works if the output file is NOT saved in   !!!"
    echo "!!!                     binary format                           !!!"
    echo "-------------------------------------------------------------------"
    echo "<APD_LINE> guide:"
    echo "  VTP/VTU --> Line 6"
    echo "  VTK     --> Last line"
    echo "-------------------------------------------------------------------"
    exit 1
fi

for i in `ls -1v ${APD_DIR}/apd-map*`; do sed -n "${APD_LINE}p" $i | awk -v "col=${CELL_ID}" -F ' ' '{print $col}'  ; done > ${APD_OUT}

