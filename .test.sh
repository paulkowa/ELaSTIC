#!/bin/bash

clear

rm -rf build/*
cd build/
cmake ../ -DCMAKE_INSTALL_PREFIX=. -DCMAKE_CXX_FLAGS="-O3"
make -j6 install
