#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task3_ts
#SBATCH -o %x.out -e %x.err
#SBATCH -N 1 -c 20
#SBATCH -t 0-00:30:00

g++ task3.cpp msort.cpp -Wall -O3 -std=c++17 -o task3 -fopenmp

for i in {1..10}; do
  ts=$((2 ** i))
  echo "Running with ts = $ts"
  ./task3 1000000 8 $ts
  echo
done
