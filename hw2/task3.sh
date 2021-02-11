#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task3
#SBATCH -o task3.out -e task3.err

if [ "$1" = "compile" ]; then
  g++ matmul.cpp task3.cpp -Wall -O3 -o task3
elif [ "$1" = "run" ]; then
  ./task3 "$2"
elif [ "$1" = "clean" ]; then
  rm -f task3 task3.out task3.err
else
  echo "./task3.sh [compile | run | clean]"
fi
