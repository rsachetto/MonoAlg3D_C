#!/bin/bash

# !!! PUT ME INTO THE ROOT FOLDER IN ORDER TO RUN ---> MonoAlg3D_C/

# CHANGE THOSE PATHS IF YOU WANT ...
PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/"
CONFIG_FILES_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/elnaz_samples/paper-purkinje-coupled"

#echo "=======================================================SCENARIO 1============================================================"
#${PROGRAM_PATH}/bin/MonoAlg3D_C -c ${CONFIG_FILES_PATH}/elnaz_purkinje_coupled_sc1.ini

#echo "=======================================================SCENARIO 2============================================================"
for NUMBER in $(seq 1 20); do
    #echo "[!] Calculating error for individual 1_${NUMBER} ..."
    ${PROGRAM_PATH}/bin/MonoAlg3D -c ${CONFIG_FILES_PATH}/elnaz_purkinje_coupled_sc2_${NUMBER}.ini
    #echo
done

#echo "=======================================================SCENARIO 3============================================================"
for NUMBER in $(seq 1 20); do
    #echo "[!] Calculating error for individual 2_${NUMBER} ..."
    ${PROGRAM_PATH}/bin/MonoAlg3D -c ${CONFIG_FILES_PATH}/elnaz_purkinje_coupled_sc3_${NUMBER}.ini
    #echo
done
