#!/bin/bash

GET_APS_SCRIPT_PATH="../../"
CELL_INDEX="6240"

#echo "[!] Working on Scenario 0 ..."
#$GET_APS_SCRIPT_PATH/getAps.sh ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC0 V 6 ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC0/cell-$CELL_INDEX.txt $CELL_INDEX
#cp ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC0/cell-$CELL_INDEX.txt aps/cell-$CELL_INDEX-sc0.txt

echo "[!] Working on Scenario 1.1 ..."
$GET_APS_SCRIPT_PATH/getAps.sh ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC1.1 V 6 ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC1.1/cell-$CELL_INDEX.txt $CELL_INDEX
cp ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC1.1/cell-$CELL_INDEX.txt aps/cell-$CELL_INDEX-sc1_1.txt

echo "[!] Working on Scenario 1.2 ..."
$GET_APS_SCRIPT_PATH/getAps.sh ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC1.2 V 6 ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC1.2/cell-$CELL_INDEX.txt $CELL_INDEX
cp ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC1.2/cell-$CELL_INDEX.txt aps/cell-$CELL_INDEX-sc1_2.txt

echo "[!] Working on Scenario 2.1 ..."
$GET_APS_SCRIPT_PATH/getAps.sh ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC2.1 V 6 ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC2.1/cell-$CELL_INDEX.txt $CELL_INDEX
cp ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC2.1/cell-$CELL_INDEX.txt aps/cell-$CELL_INDEX-sc2_1.txt

echo "[!] Working on Scenario 2.2 ..."
$GET_APS_SCRIPT_PATH/getAps.sh ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC2.2 V 6 ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC2.2/cell-$CELL_INDEX.txt $CELL_INDEX
cp ../../../outputs/elnaz_plain_mixed_models_1cm_original_steadystate_SC2.2/cell-$CELL_INDEX.txt aps/cell-$CELL_INDEX-sc2_2.txt
