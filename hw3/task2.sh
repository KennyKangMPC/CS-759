#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task2
#SBATCH -o task2.out -e task2.err
#SBATCH --gres=gpu:1

if [ "$1" = "compile" ]; then
  module load cuda
  nvcc task2.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -o task2
elif [ "$1" = "run" ]; then
  ./task2
elif [ "$1" = "clean" ]; then
  rm -f task2 task2.out task2.err
else
  echo "./$0 [compile | run | clean]"
fi
