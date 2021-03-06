#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task1-512
#SBATCH -o task1-512.out -e task1-512.err
#SBATCH --gres=gpu:1
#SBATCH -t 0-00:30:00

module load cuda
nvcc task1.cu reduce.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -o task1
for i in {10..30}; do
  n=$((2 ** i))
  echo "Running with n = 2^$i = $n"
  ./task1 $n 512
  echo
done