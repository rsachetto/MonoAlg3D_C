#!/bin/bash

./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_cco_seed-1562005512_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu
./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_cco_seed-1562009134_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu
./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_cco_seed-1562013988_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu
./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_cco_seed-1562042299_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu
./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_cco_seed-1562046115_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu

echo
./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_co_seed-1562005512_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu
./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_co_seed-1562009134_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu
./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_co_seed-1562013988_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu
./bin/ErrorCalculator ../../outputs/elizabeth_coupled_arpf_gold_LV/activation-map.vtu ../../outputs/elizabeth_coupled_arpf_co_seed-1562046115_offset-9_LV/activation-map.vtu elizabeth_cco_error_at.vtu


