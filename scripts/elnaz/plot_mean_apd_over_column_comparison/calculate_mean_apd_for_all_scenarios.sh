#!/bin/bash

echo "[!] Working on Scenario 0 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc0.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc0.txt

echo "[!] Working on Scenario 1 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc1.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc1.txt

echo "[!] Working on Scenario 2 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc2.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc2.txt

echo "[!] Working on Scenario 2.1 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc2_1.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc2_1.txt

echo "[!] Working on Scenario 2.2 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc2_2.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc2_2.txt

echo "[!] Working on Scenario 3 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc3.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc3.txt

echo "[!] Working on Scenario 3.1 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc3_1.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc3_1.txt

echo "[!] Working on Scenario 3.2 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc3_2.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc3_2.txt

echo "[!] Working on Scenario 3.3 ..."
./bin/ApdPlot inputs/apd-map-steadystate-sc3_3.vtu 10000 100
mv ./outputs/mean_apd.txt ./outputs/mean_apd_sc3_3.txt