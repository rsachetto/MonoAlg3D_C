#!/bin/bash

echo "[!] Working on Scenario 0 ..."
#./bin/ApdPlot inputs/apd-map-steadystate-sc0.vtu 10000 100
./bin/ApdPlot inputs/apd-map-cable-sc0.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc0.txt

#echo "[!] Working on Scenario 1.1 ..."
#./bin/ApdPlot inputs/apd-map-steadystate-sc1_1.vtu 10000 100
#mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc1_1.txt

#echo "[!] Working on Scenario 1.2 ..."
#./bin/ApdPlot inputs/apd-map-steadystate-sc1_2.vtu 10000 100
#mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc1_2.txt

#echo "[!] Working on Scenario 2.1 ..."
#./bin/ApdPlot inputs/apd-map-steadystate-sc2_1.vtu 10000 100
#mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc2_1.txt

#echo "[!] Working on Scenario 2.2 ..."
#./bin/ApdPlot inputs/apd-map-steadystate-sc2_2.vtu 10000 100
#mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc2_2.txt

