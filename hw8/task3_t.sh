#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task3_t
#SBATCH -o %x.out -e %x.err
#SBATCH -N 1 -c 20
#SBATCH -t 0-00:30:00

g++ task3.cpp msort.cpp -Wall -O3 -std=c++17 -o task3 -fopenmp

for t in {1..20}; do
  echo "Running with t = $t"
  ./task3 1000000 $t 32
  echo
done
