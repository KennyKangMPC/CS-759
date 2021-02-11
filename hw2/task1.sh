#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task1
#SBATCH -o task1.out -e task1.err

if [ "$1" = "compile" ]; then
  g++ scan.cpp task1.cpp -Wall -O3 -o task1
elif [ "$1" = "run" ]; then
  for i in {10..30}; do
    n=$((2 ** i))
    echo "Running with n = 2^$i = $n"
    ./task1 $n
    echo
  done
elif [ "$1" = "plot" ]; then
  python task1-plot.py
elif [ "$1" = "clean" ]; then
  rm -f task1 task1.err task1.out task1.pdf
else
  echo "./task1.sh [compile | run | plot | clean]"
fi
