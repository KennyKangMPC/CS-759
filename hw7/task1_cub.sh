#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task1_cub
#SBATCH -o %x.out -e %x.err
#SBATCH --gres=gpu:1
#SBATCH -t 0-00:30:00

module load cuda
nvcc task1_cub.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -std c++17 -o task1_cub

for i in {10..30}; do
  n=$((2 ** i))
  echo "Running with n = 2^$i = $n"
  ./task1_cub $n
  echo
done
