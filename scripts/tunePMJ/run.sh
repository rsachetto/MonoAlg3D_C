#!/bin/bash

# =====================================================================================================================
# This script is responsible for tuning the PMJ delay parameters
# =====================================================================================================================

# GLOBALS
PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/scripts/tunePMJ"
CONFIG_FILE_PATH="/home/berg/Github/MonoAlg3D_C/scripts/tunePMJ/configs"
OUTPUT_FILE_PATH="/home/berg/Github/MonoAlg3D_C/scripts/tunePMJ/outputs/purkinje_cuboid"
RPMJS=( 2.1 2.2 )

rm pmj_delay.txt
rm rpmj.txt

for RPMJ in "${RPMJS[@]}"; do
    ${PROGRAM_PATH}/bin/tunePMJ ${RPMJ}
    echo ${RPMJ} >> rpmj.txt
    grep "\[purkinje_coupling\]" ${OUTPUT_FILE_PATH}/outputlog.txt | cut -d' ' -f24 >> pmj_delay.txt
done

paste rpmj.txt pmj_delay.txt > rpmj_delay.txt
