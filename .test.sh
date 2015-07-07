#!/bin/bash

clear

rm -rf build/*
cd build/
cmake ../ -DCMAKE_INSTALL_PREFIX=.
make -j6 install
