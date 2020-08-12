#!/bin/bash

SEEDS=( 1562002891 1562002894 1562005513 1562005553 1562006177 1562007596 1562008172 1562008424 1562009135 1562009769 )

# ========================================================================================================================
# 0) REFERENCE
# ========================================================================================================================
#MONOALG_PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D"
#MONOALG_INPUT_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/05_Lucas/02_Oxford/00_Reference"

#${MONOALG_PROGRAM_PATH} -c ${MONOALG_INPUT_PATH}/simple_reference.ini

# ========================================================================================================================
# 1) MINIMIZE TERMINALS
# ========================================================================================================================
#MONOALG_PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D"
#MONOALG_INPUT_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/05_Lucas/02_Oxford/01_Minimize_Terminals"

#for SEED in "${SEEDS[@]}"; do
#    ${MONOALG_PROGRAM_PATH} -c ${MONOALG_INPUT_PATH}/min:terminals_seed-${SEED}.ini
#done

# ========================================================================================================================
# 2) MINIMIZE TOTAL
# ========================================================================================================================
#MONOALG_PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D"
#MONOALG_INPUT_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/05_Lucas/02_Oxford/02_Minimize_Total"

#for SEED in "${SEEDS[@]}"; do
#    ${MONOALG_PROGRAM_PATH} -c ${MONOALG_INPUT_PATH}/min:total_seed-${SEED}.ini
#done

# ========================================================================================================================
# 3) MINIMIZE TERMINALS LINKED
# ========================================================================================================================
#MONOALG_PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D"
#MONOALG_INPUT_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/05_Lucas/02_Oxford/03_Minimize_Terminals_Linked"

#for SEED in "${SEEDS[@]}"; do
#    ${MONOALG_PROGRAM_PATH} -c ${MONOALG_INPUT_PATH}/min:terminals-linked_seed-${SEED}.ini
#done

# ========================================================================================================================
# 4) MINIMIZE TOTAL LINKED
# ========================================================================================================================
#MONOALG_PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D"
#MONOALG_INPUT_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/05_Lucas/02_Oxford/04_Minimize_Total_Linked"

#for SEED in "${SEEDS[@]}"; do
#    ${MONOALG_PROGRAM_PATH} -c ${MONOALG_INPUT_PATH}/min:total-linked_seed-${SEED}.ini
#done

# ========================================================================================================================
# 5) MINIMIZE TERMINALS LINKED SOFT
# ========================================================================================================================
#MONOALG_PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D"
#MONOALG_INPUT_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/05_Lucas/02_Oxford/05_Minimize_Terminals_Linked_Soft_Pruning"

#for SEED in "${SEEDS[@]}"; do
#    ${MONOALG_PROGRAM_PATH} -c ${MONOALG_INPUT_PATH}/min:terminals-linked-soft_seed-${SEED}.ini
#done

# ========================================================================================================================
# 6) MINIMIZE TERMINALS LINKED MODERATE
# ========================================================================================================================
MONOALG_PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D"
MONOALG_INPUT_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/05_Lucas/02_Oxford/06_Minimize_Terminals_Linked_Moderate_Pruning"

for SEED in "${SEEDS[@]}"; do
    ${MONOALG_PROGRAM_PATH} -c ${MONOALG_INPUT_PATH}/min:terminals-linked-moderate_seed-${SEED}.ini
done

# ========================================================================================================================
# 7) MINIMIZE TERMINALS LINKED HEAVY
# ========================================================================================================================
MONOALG_PROGRAM_PATH="/home/berg/Github/MonoAlg3D_C/bin/MonoAlg3D"
MONOALG_INPUT_PATH="/home/berg/Github/MonoAlg3D_C/example_configs/05_Lucas/02_Oxford/07_Minimize_Terminals_Linked_Heavy_Pruning"

for SEED in "${SEEDS[@]}"; do
    ${MONOALG_PROGRAM_PATH} -c ${MONOALG_INPUT_PATH}/min:terminals-linked-heavy_seed-${SEED}.ini
done

