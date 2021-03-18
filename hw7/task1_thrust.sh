#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task1_thrust
#SBATCH -o %x.out -e %x.err
#SBATCH --gres=gpu:1
#SBATCH -t 0-00:30:00

module load cuda
nvcc task1_thrust.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -o task1_thrust

for i in {10..30}; do
  n=$((2 ** i))
  echo "Running with n = 2^$i = $n"
  ./task1_thrust $n
  echo
done
