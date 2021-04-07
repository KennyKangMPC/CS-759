#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task2_pure_omp
#SBATCH -o %x.out -e %x.err
#SBATCH --nodes=1 --cpus-per-task=20 --ntasks-per-node=1
#SBATCH -t 0-00:30:00

module load mpi/openmpi
module load gcc/10.2.0

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

g++ task2_pure_omp.cpp reduce.cpp -Wall -O3 -o task2_pure_omp -fopenmp -fno-tree-vectorize -march=native -fopt-info-vec

n=20000000
for t in {1..20}; do
  echo "Running with n = $n, t = $t"
  ./task2_pure_omp $n $t
  echo
done

