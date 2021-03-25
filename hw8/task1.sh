#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task1
#SBATCH -o %x.out -e %x.err
#SBATCH -N 1 -c 20
#SBATCH -t 0-00:30:00

g++ task1.cpp matmul.cpp -Wall -O3 -o task1 -fopenmp

for t in {1..20};
do
  echo "Running with t = $t"
  ./task1 1024 $t
  echo
done
