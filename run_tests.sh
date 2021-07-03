#!/bin/bash

INITIAL_PARAMS=""
VALID_TESTS="libs mesh solvers simulations"

function PRINT_USAGE() {
    echo "Usage $0 [test]" >&2;
    echo "Valid tests: all ${VALID_TESTS} (default is all)" >&2;
}

function RUN_TEST() {
    lower_case_name=$1
    test_bin_name="$(tr '[:lower:]' '[:upper:]' <<< "${lower_case_name:0:1}")${lower_case_name:1}"
    test_bin_name="./tests_bin/Test${test_bin_name}"

    if [ "$lower_case_name" == "all" ]; then
        for t in ${VALID_TESTS}; do
            t_bin_name="$(tr '[:lower:]' '[:upper:]' <<< "${t:0:1}")${t:1}"
            echo "RUNNING ./tests_bin/Test${t_bin_name}"
            ./tests_bin/Test"${t_bin_name}" --jobs 1
        done
    else
        if [ -f "${test_bin_name}" ]; then
            "${test_bin_name}" --jobs 1
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

RUN_TEST "$TESTS"
