#!/usr/bin/env bash
AP_DIR=$1
AP_PREFIX=$2
AP_LINE=$3
PRINT_RATE=$4

if [[ "$#" -ne 4 ]]; then
    echo "-------------------------------------------------------------------"
    echo "Usage:> $0 <AP_DIR> <AP_PREFIX> <AP_LINE> <PRINT_RATE>"
    echo "-------------------------------------------------------------------"
    echo "<AP_DIR> = Directory where the results of the simulation are stored"
    echo "<AP_PREFIX> = Prefix of the results files"
    echo "<AP_LINE> = Line where the action potential information is stored"
    echo "<PRINT_RATE> = The print rate of the simulation"
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

#for i in `ls -1v ${AP_DIR}/${AP_PREFIX}*`; do sed -n "${AP_LINE}p" $i | awk -v "col=${CELL_ID}" -F ' ' '{print $col}'  ; done > ${AP_OUT}

TIMESTEP=0
for i in `ls -1v ${AP_DIR}/${AP_PREFIX}*`; do
    AP_OUT="aps/timestep-$TIMESTEP.txt"
    sed -n "${AP_LINE}p" $i | awk -F ' ' '{ for (i=0; i < 10000; i++) { printf "%g\n", $i } }' > ${AP_OUT}
    let "TIMESTEP=TIMESTEP+$PRINT_RATE"
done

let "TOTAL_TIMESTEPS=TIMESTEP-$PRINT_RATE"

CMD=""
for i in $(seq 0 $PRINT_RATE $TOTAL_TIMESTEPS); do
    CMD="$CMD aps/timestep-$i.txt"
done
#echo $CMD
paste $CMD > timesteps.txt
