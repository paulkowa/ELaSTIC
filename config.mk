# Compiler
CXX=g++
AR=ar

# MPI Compiler
MPICXX=mpicxx

# OpenMP switch
OMPFLAGS=-fopenmp

# Basic compiler flags
CXXFLAGS=-O3

# Boost library (if not in the standard location)
#BOOST_INCLUDE=-I/opt/boost/include
#BOOST_LIB=-L/opt/boost/lib

# Use MPE (for profiling only!)
# MPE library (if needed, and not in standard location)
#MPE_LIB=-L/opt/mpe2/lib
#MPE_FLAGS=-llmpe -lmpe
#MPE_INCLUDE=-I/opt/mpe2/include
#WITH_MPE=-DWITH_MPE
