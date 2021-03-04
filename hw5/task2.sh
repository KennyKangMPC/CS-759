#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task2
#SBATCH -o task2.out -e task2.err
#SBATCH --gres=gpu:1
#SBATCH -t 0-00:30:00

module load cuda
nvcc task2.cu matmul.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -o task2
for i in {5..15}; do
  n=$((2 ** i))
  echo "Running with n = 2^$i = $n"
  ./task2 $n 32
  echo
done
