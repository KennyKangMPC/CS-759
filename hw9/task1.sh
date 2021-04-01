#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task1
#SBATCH -o %x.out -e %x.err
#SBATCH -N 1 -c 20
#SBATCH -t 0-00:30:00

g++ task1.cpp cluster.cpp -Wall -O3 -std=c++17 -o task1 -fopenmp

for t in {1..10}; do
  echo "Running with t = $t"
  ./task1 5040000 $t
  echo
done
