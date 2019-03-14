#!/bin/bash

BUILD_DIR="build"
BUILD_TYPE="Release"

if [[ "$#" -eq 1 ]]; then
    BUILD_DIR=$1
fi

if [[ "$#" -eq 2 ]]; then
    BUILD_DIR=$1
    BUILD_TYPE=$2
fi


if [[ ! -d "${BUILD_DIR}" ]]; then
  echo "Directory ${BUILD_DIR} does not exist. Creating."
  mkdir ${BUILD_DIR}
fi

cd ${BUILD_DIR}; cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ..; make
