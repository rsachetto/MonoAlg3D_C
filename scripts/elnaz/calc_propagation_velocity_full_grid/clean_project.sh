#!/bin/bash

BUILD_DIR="build"

if [[ "$#" -eq 1 ]]; then
    BUILD_DIR=$1
fi

if [[ ! -d "${BUILD_DIR}" ]]; then
  echo "Directory ${BUILD_DIR} does not exist"
  exit
fi

rm -fr ${BUILD_DIR}/*
rm -fr bin/*
