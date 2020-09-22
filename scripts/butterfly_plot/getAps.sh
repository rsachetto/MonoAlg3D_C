#!/usr/bin/env bash
AP_DIR=$1
AP_PREFIX=$2
AP_LINE=$3
TOTAL_NUMBER_CELLS=$4
PRINT_RATE=$5
SIMULATION_NUMBER=$6

if [[ "$#" -ne 6 ]]; then
    echo "---------------------------------------------------------------------------------"
    echo "Usage:> $0 <AP_DIR> <AP_PREFIX> <AP_LINE> <TOTAL_NUMBER_CELLS> <PRINT_RATE> <SIMULATION_NUMBER>"
    echo "---------------------------------------------------------------------------------"
    echo "<AP_DIR> = Directory where the results of the simulation are stored"
    echo "<AP_PREFIX> = Prefix of the results files"
    echo "<AP_LINE> = Line where the action potential information is stored"
    echo "<TOTAL_NUMBER_CELLS> = Total number of grid cells"
    echo "<PRINT_RATE> = The print rate of the simulation"
    echo "<SIMULATION_NUMBER> = Simulation number"
    echo "---------------------------------------------------------------------------------"
    echo "!!! This script only works if the output file is NOT saved in   !!!"
    echo "!!!                     binary format                           !!!"
    echo "---------------------------------------------------------------------------------"
    echo "<AP_LINE> guide:"
    echo "  VTP/VTU --> Line 6"
    echo "  VTK     --> Last line"
    echo "---------------------------------------------------------------------------------"
    exit 1
fi

# Get the transmembrane potential from each of the .vtu files and dump it into a .txt
TIMESTEP=0
for i in `ls -1v ${AP_DIR}/${AP_PREFIX}*`; do
    AP_OUT="aps/timestep-$SIMULATION_NUMBER-$TIMESTEP.txt"
    sed -n "${AP_LINE}p" $i | awk -v TOTAL_NUMBER_CELLS="$TOTAL_NUMBER_CELLS" -F ' ' '{ for (i=0; i < TOTAL_NUMBER_CELLS; i++) { printf "%g\n", $i } }' > ${AP_OUT}
    let "TIMESTEP=TIMESTEP+$PRINT_RATE"
done

# Adjust the total number of timesteps
let "TOTAL_TIMESTEPS=TIMESTEP-$PRINT_RATE"

# Concatenate the filename from each timestep
CMD=""
for i in $(seq 0 $PRINT_RATE $TOTAL_TIMESTEPS); do
    CMD="$CMD aps/timestep-$SIMULATION_NUMBER-$i.txt"
done

# Use the 'paste' command to append each timestep into a column on the output file
paste $CMD > timesteps-$SIMULATION_NUMBER.txt
