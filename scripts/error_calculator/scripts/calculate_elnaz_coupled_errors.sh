#!/bin/bash

# CHANGE THOSE PATHS IF YOU WANT ...
PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/scripts/error_calculator"
OUTPUTS_PATH="/home/berg/Github/MonoAlg3D_C/outputs"
ERROR_OUTPUT_PATH="/home/berg/Github/MonoAlg3D_C/scripts/error_calculator/outputs/elnaz-purkinje-coupled-errors"

#echo "=======================================================SCENARIO 1============================================================"
for NUMBER in $(seq 1 20); do
    #echo "[!] Calculating error for individual 1_${NUMBER} ..."
    ${PROGRAM_PATH}/bin/ErrorCalculator ${OUTPUTS_PATH}/elnaz_purkinje_coupled_sc0/activation-map.vtu ${OUTPUTS_PATH}/elnaz_purkinje_coupled_sc1_${NUMBER}/activation-map.vtu ${ERROR_OUTPUT_PATH}/error-activation-map-sc0-sc1_${NUMBER}.vtu >> ../outputs/results-error-activation-map-sc1.txt
    #echo
done

#echo "=======================================================SCENARIO 2============================================================"
for NUMBER in $(seq 1 20); do
    #echo "[!] Calculating error for individual 2_${NUMBER} ..."
    ${PROGRAM_PATH}/bin/ErrorCalculator ${OUTPUTS_PATH}/elnaz_purkinje_coupled_sc0/activation-map.vtu ${OUTPUTS_PATH}/elnaz_purkinje_coupled_sc2_${NUMBER}/activation-map.vtu ${ERROR_OUTPUT_PATH}/error-activation-map-sc0-sc2_${NUMBER}.vtu >> ../outputs/results-error-activation-map-sc2.txt
    #echo
done

#echo "=======================================================SCENARIO 1============================================================"
for NUMBER in $(seq 1 20); do
    #echo "[!] Calculating error for individual 1_${NUMBER} ..."
    ${PROGRAM_PATH}/bin/ErrorCalculator ${OUTPUTS_PATH}/elnaz_purkinje_coupled_sc0/apd-map.vtu ${OUTPUTS_PATH}/elnaz_purkinje_coupled_sc1_${NUMBER}/apd-map.vtu ${ERROR_OUTPUT_PATH}/error-apd-map-sc0-sc1_${NUMBER}.vtu >> ../outputs/results-error-apd-map-sc1.txt
    #echo
done

#echo "=======================================================SCENARIO 2============================================================"
for NUMBER in $(seq 1 20); do
    #echo "[!] Calculating error for individual 1_${NUMBER} ..."
    ${PROGRAM_PATH}/bin/ErrorCalculator ${OUTPUTS_PATH}/elnaz_purkinje_coupled_sc0/apd-map.vtu ${OUTPUTS_PATH}/elnaz_purkinje_coupled_sc2_${NUMBER}/apd-map.vtu ${ERROR_OUTPUT_PATH}/error-apd-map-sc0-sc2_${NUMBER}.vtu >> ../outputs/results-error-apd-map-sc2.txt
    #echo
done
