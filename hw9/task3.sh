#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task3
#SBATCH -o %x.out -e %x.err
#SBATCH --ntasks-per-node=2
#SBATCH -t 0-00:30:00

module load mpi/openmpi

mpicxx task3.cpp -Wall -O3 -o task3

for i in {1..25}; do
  n=$((2 ** i))
  echo "Running with n = 2^$i = $n"
  mpirun -np 2 task3 $n
  echo
done
