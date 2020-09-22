#!/bin/bash

for i in $(seq 1 3); do 
	./bin/ErrorCalculator inputs/corner-activation-map-sc0.vtu inputs/corner-activation-map-sc$i.vtu outputs/corner-error-activation-time-sc0-sc$i.vtu
	./bin/ErrorCalculator inputs/corner-conduction-velocity-sc0.vtu inputs/corner-conduction-velocity-sc$i.vtu outputs/corner-error-conduction-velocity-sc0-sc$i.vtu
done


