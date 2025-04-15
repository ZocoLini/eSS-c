#!/bin/bash

source "$HOME/intel/oneapi/setvars.sh"

rm -rf build
mkdir build

(
  cd build || exit 1
  cmake ..
  make
)