#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task2
#SBATCH -o %x.out -e %x.err
#SBATCH -N 1 -c 20
#SBATCH -t 0-00:30:00

g++ task2.cpp convolution.cpp -Wall -O3 -std=c++17 -o task2 -fopenmp

for t in {1..20}; do
  echo "Running with t = $t"
  ./task2 1024 $t
  echo
done
