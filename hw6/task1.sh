#!/usr/bin/env bash

#SBATCH -p ppc --gres=gpu:v100:1 -t 0-00:20:00
#SBATCH -J task1
#SBATCH -o task1.out -e task1.err

#module purge
#module load cuda/10.1 clang/7.0.0
#nvcc task1.cu mmul.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -lcublas -ccbin $CC -o task1
for i in {5..15}; do
  n=$((2 ** i))
  echo "Running with n = 2^$i = $n"
  ./task1 $n 6
  echo
done
