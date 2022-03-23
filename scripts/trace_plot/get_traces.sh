#!/bin/bash
# Author: Lucas Berg
# This script is used to get the traces for cable simulation related to the:
# "Performance Improvement of the ToRORd model"  and
# "Performance Improvement of the Trovato model"

STATE_VECTOR_EULER_SACHETTO_FILEPATH="/home/berg/Github/MonoAlg3D_C/outputs/trovato_cable_EulerAdapt_Sachetto_cpu"
STATE_VECTOR_EULER_JHONNY_FILEPATH="/home/berg/Github/MonoAlg3D_C/outputs/trovato_cable_EulerAdapt_Jhonny_cpu"
STATE_VECTOR_RL_JHONNY_FILEPATH="/home/berg/Github/MonoAlg3D_C/outputs/trovato_cable_RLAdapt_Jhonny_cpu"
# ToRORd-dynCl
#STATE_VECTOR_NAMES=( "v" "CaMKt" "cass" "nai" "nass" "ki" "kss" "cansr" "cajsr" "cai" "m" "h" "j" "hp" "jp" "mL" "hL" "hLp" "a" "iF" "iS" "ap" "iFp" "iSp" "d" "ff" "fs" "fcaf" "fcas" "jca" "ffp" "fcafp" "nca_ss" "nca_i" "C1" "C2" "C3" "I" "O" "xs1" "xs2" "Jrel_np" "Jrel_p" )
# Trovato
STATE_VECTOR_NAMES=( "v" "CaMKt" "cass" "nai" "nasl" "nass" "ki" "kss" "ksl" "cai" "casl" "cansr" "cajsr" "cacsr" "Jrel1" "Jrel2" "m" "hf" "hs" "j" "hsp" "jp" "mL" "hL" "hLp" "a" "i1" "i2" "d" "ff" "fs" "fcaf" "fcas" "jca" "ffp" "fcafp" "nca" "b" "g" "xrf" "xs1" "xs2" "y" )

#COUNTER=2
#for NAME in "${STATE_VECTOR_NAMES[@]}"; do
#    cut -f${COUNTER} -d' ' ${STATE_VECTOR_EULER_SACHETTO_FILEPATH}/sv_0.txt > ${STATE_VECTOR_EULER_SACHETTO_FILEPATH}/${NAME}.txt
#    let COUNTER=${COUNTER}+2
#done

#COUNTER=2
#for NAME in "${STATE_VECTOR_NAMES[@]}"; do
#    cut -f${COUNTER} -d' ' ${STATE_VECTOR_EULER_JHONNY_FILEPATH}/sv_0.txt > ${STATE_VECTOR_EULER_JHONNY_FILEPATH}/${NAME}.txt
#    let COUNTER=${COUNTER}+2
#done

#COUNTER=2
#for NAME in "${STATE_VECTOR_NAMES[@]}"; do
#    cut -f${COUNTER} -d' ' ${STATE_VECTOR_RL_JHONNY_FILEPATH}/sv_0.txt > ${STATE_VECTOR_RL_JHONNY_FILEPATH}/${NAME}.txt
#    let COUNTER=${COUNTER}+2
#done

for NAME in "${STATE_VECTOR_NAMES[@]}"; do
    python plot_comparison_sv.py ${STATE_VECTOR_EULER_SACHETTO_FILEPATH}/${NAME}.txt ${STATE_VECTOR_EULER_JHONNY_FILEPATH}/${NAME}.txt ${STATE_VECTOR_RL_JHONNY_FILEPATH}/${NAME}.txt 0.02 1 ${NAME}
done


