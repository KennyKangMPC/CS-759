#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task1
#SBATCH -o task1.out -e task1.err
#SBATCH --gres=gpu:1

if [ "$1" = "compile" ]; then
  module load cuda
  nvcc task1.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -o task1
elif [ "$1" = "run" ]; then
  ./task1
elif [ "$1" = "clean" ]; then
  rm -f task1 task1.out task1.err
else
  echo "./$0 [compile | run | clean]"
fi
