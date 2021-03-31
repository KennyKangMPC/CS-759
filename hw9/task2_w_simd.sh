#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task2_w_simd
#SBATCH -o %x.out -e %x.err
#SBATCH -N 1 -c 20
#SBATCH -t 0-00:30:00

module load gcc/10.2.0

g++ -DSIMD task2.cpp montecarlo.cpp -Wall -O3 -std=c++17 -o task2 -fopenmp -fno-tree-vectorize -march=native -fopt-info-vec

for t in {1..10}; do
  echo "Running with t = $t"
  ./task2 1000000 $t
  echo
done
