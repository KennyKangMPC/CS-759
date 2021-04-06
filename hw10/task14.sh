#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task14
#SBATCH -o %x.out -e %x.err
#SBATCH -N 1 -c 1
#SBATCH -t 0-00:30:00

MACRO_FLAGS="-DDATA_T=float -DOP=* -DIDENT=1.f"
N=1000000

module load gcc/10.2.0

g++ $MACRO_FLAGS task1.cpp optimize.cpp -Wall -O3 -std=c++17 -o task1 -fno-tree-vectorize -ffast-math
./task1 $N

g++ $MACRO_FLAGS task1.cpp optimize.cpp -Wall -O3 -std=c++17 -o task1 -march=native -fopt-info-vec -ffast-math
./task1 $N
