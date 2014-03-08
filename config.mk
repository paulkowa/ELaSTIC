# -----------------------------------------------------

# Prefix of the directory where ELaSTIC should be installed
PREFIX=/usr/local

# Change your preferred C++ compiler
CXX=g++

# Change your preferred MPI compiler
# Most likely no changes needed
MPICXX=mpicxx

# Change OpenMP switch for your C++ compiler
# Typical values (http://openmp.org/wp/openmp-compilers/):
# -fopenmp for GCC (g++) and Clang/LLVM (clang++)
# -openmp for Intel (icpc)
# -mp for Portland Group (pgCC)
OMPFLAGS=-fopenmp

# Tune C++ compiler optimization flags
# Most likely no changes needed
# Remove -std=c++0x for older compilers
CXXFLAGS=-std=c++0x -O3

# Set the Boost C++ Libraries paths if different from standard
#BOOST_INCLUDE=-I/opt/boost/include
#BOOST_LIB=-L/opt/boost/lib

# Most likely AR does not have to be changed
AR=ar

# -----------------------------------------------------

# This section is for developers
# Enable MPE for profiling only

# MPE library if required, and not in standard location
#MPE_LIB=-L/opt/mpe2/lib
#MPE_FLAGS=-llmpe -lmpe
#MPE_INCLUDE=-I/opt/mpe2/include
#WITH_MPE=-DWITH_MPE
