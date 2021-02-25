#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task1
#SBATCH -o task1.out -e task1.err
#SBATCH --gres=gpu:1
#SBATCH -t 0-00:30:00

nvcc task1.cu matmul.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -o task1
for i in {5..14}; do
  n=$((2 ** i))
  echo "Running with n = 2^$i = $n"
  ./task1 $n 1024
  echo
done
