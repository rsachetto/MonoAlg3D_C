#!/bin/bash

# CHANGE THOSE PATHS IF YOU WANT ...
OUTPUT_DIR_PATH="/home/berg/Github/MonoAlg3D_C/outputs"
DATA_FOLDER_PATH="/home/berg/Github/MonoAlg3D_C/scripts/elnaz/scientific-reports-purkinje"
EXTRACT_CELL_DATA_PROGRAM="/home/berg/Github/MonoAlg3D_C/scripts/elnaz/extract-celldata"

#=======================================================SCENARIO 1============================================================
cp ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc1/activation-map.vtu ${DATA_FOLDER_PATH}/sc1/activation-map-sc1.vtu
${EXTRACT_CELL_DATA_PROGRAM}/bin/ExtractCellData ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc1/activation-map.vtu ${DATA_FOLDER_PATH}/sc1/activation-map-sc1.txt

cp ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc1/apd-map.vtu ${DATA_FOLDER_PATH}/sc1/apd-map-sc1.vtu
${EXTRACT_CELL_DATA_PROGRAM}/bin/ExtractCellData ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc1/apd-map.vtu ${DATA_FOLDER_PATH}/sc1/apd-map-sc1.txt


#=======================================================SCENARIO 2============================================================
for NUMBER in $(seq 1 20); do
    cp ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc2_${NUMBER}/activation-map.vtu ${DATA_FOLDER_PATH}/sc2/activation-map-sc2_${NUMBER}.vtu
    ${EXTRACT_CELL_DATA_PROGRAM}/bin/ExtractCellData ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc2_${NUMBER}/activation-map.vtu ${DATA_FOLDER_PATH}/sc2/activation-map-sc2_${NUMBER}.txt

    cp ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc2_${NUMBER}/apd-map.vtu ${DATA_FOLDER_PATH}/sc2/apd-map-sc2_${NUMBER}.vtu
    ${EXTRACT_CELL_DATA_PROGRAM}/bin/ExtractCellData ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc2_${NUMBER}/apd-map.vtu ${DATA_FOLDER_PATH}/sc2/apd-map-sc2_${NUMBER}.txt
done

#=======================================================SCENARIO 2============================================================
for NUMBER in $(seq 1 20); do
    cp ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc3_${NUMBER}/activation-map.vtu ${DATA_FOLDER_PATH}/sc3/activation-map-sc3_${NUMBER}.vtu
    ${EXTRACT_CELL_DATA_PROGRAM}/bin/ExtractCellData ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc3_${NUMBER}/activation-map.vtu ${DATA_FOLDER_PATH}/sc3/activation-map-sc3_${NUMBER}.txt

    cp ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc3_${NUMBER}/apd-map.vtu ${DATA_FOLDER_PATH}/sc3/apd-map-sc3_${NUMBER}.vtu
    ${EXTRACT_CELL_DATA_PROGRAM}/bin/ExtractCellData ${OUTPUT_DIR_PATH}/elnaz_purkinje_coupled_sc3_${NUMBER}/apd-map.vtu ${DATA_FOLDER_PATH}/sc3/apd-map-sc3_${NUMBER}.txt
done
