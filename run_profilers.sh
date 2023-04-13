#!/bin/bash

INITIAL_PARAMS=""
VALID_TESTS="mesh custom_mesh gpusolvers solvers simulation vtu txt bin en alg"

function PRINT_USAGE() {
    echo "Usage $0 [profiler]" >&2;
    echo "Valid profilers: all ${VALID_TESTS} (default is all)" >&2;
}

function RUN_PROFILER() {
    lower_case_name=$1
    test_bin_name="$(tr '[:lower:]' '[:upper:]' <<< "${lower_case_name:0:1}")${lower_case_name:1}"
    test_bin_name="./tests_bin/${test_bin_name}Profiler"

    hardware_key=$(inxi -CNM |grep -e Topology  -e Device-1 -e Mobo | sha256sum | awk '{print $1}')

    if [ "$lower_case_name" == "all" ]; then
        for t in ${VALID_TESTS}; do
            t_bin_name="$(tr '[:lower:]' '[:upper:]' <<< "${t:0:1}")${t:1}"
            echo "-----------------------------"
            echo "Running ${t_bin_name}Profiler"
            echo "-----------------------------"
            ./tests_bin/"${t_bin_name}"Profiler "$hardware_key"
        done
    else
        if [ -f "${test_bin_name}" ]; then
            "${test_bin_name}" "$hardware_key"
        else
            PRINT_USAGE "$INITIAL_PARAMS"
        fi
    fi
}

TESTS=""

if [ "$#" -eq 0 ]; then
    TESTS="all"
elif [ "$#" -eq 1 ]; then
    TESTS="$1"
else
    PRINT_USAGE "$@"
    exit 1
fi

INITIAL_PARAMS="$*"

RUN_PROFILER "$TESTS"


