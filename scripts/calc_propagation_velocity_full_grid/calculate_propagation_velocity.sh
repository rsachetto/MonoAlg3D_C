#!/bin/bash

OUTPUT_DIR="../../outputs/plain_100_100_100_tentusscher_epi"
DT=0.02
PRINT_RATE=50
SIDE_LENGTH_X=10000
SIDE_LENGTH_Y=10000
DX=100
DY=100
ACTIVATION_MAP_PATH="inputs/corner-activation-map-sc0.vtu"

for i in $(seq 0 3); do
	ACTIVATION_MAP_PATH="inputs/corner-activation-map-sc$i.vtu"
	./bin/CalcPropagationVelocity $OUTPUT_DIR $DT $PRINT_RATE $SIDE_LENGTH_X $SIDE_LENGTH_Y $DX $DY $ACTIVATION_MAP_PATH
	mv outputs/activation_time_map.vtu outputs/corner-activation-map-sc$i.vtu
	mv outputs/conduction_velocity_map.vtu outputs/corner-conduction-velocity-sc$i.vtu
done


