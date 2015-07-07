#!/bin/bash

clear

rm -rf build/*
cd build/
cmake ../ -DCMAKE_INSTALL_PREFIX=. -DWITH_MPE=1 -DMPE_INCLUDE_DIR=/opt/mpe2-1.3.0-gnu/include -DMPE_LIB_DIR=/opt/mpe2-1.3.0-gnu/lib
make -j6 install
