#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task3
#SBATCH -o %x.out -e %x.err
#SBATCH -N 1 -c 20
#SBATCH -t 0-00:30:00

g++ task3.cpp -Wall -O3 -std=c++17 -o task3 -fopenmp
./task3
