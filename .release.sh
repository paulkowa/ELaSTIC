#!/bin/bash

rm -rf build/*
cd build/
cmake ../
make -j6
src/mpi/elastic-sketch > ../doc/elastic-sketch.txt
cd ..
rm -rf build/*
DIR=../ELaSTIC-`cat VERSION`
mkdir -p $DIR
cp -R * $DIR/
