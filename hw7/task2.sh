#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task2
#SBATCH -o %x.out -e %x.err
#SBATCH --gres=gpu:1
#SBATCH -t 0-00:30:00

module load cuda
nvcc task2.cu count.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -std c++17 -o task2

for i in {5..24}; do
  n=$((2 ** i))
  echo "Running with n = 2^$i = $n"
  ./task2 $n
  echo
done
