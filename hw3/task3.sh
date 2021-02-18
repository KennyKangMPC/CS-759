#!/usr/bin/env bash

#SBATCH -p wacc
#SBATCH -J task3
#SBATCH -o task3.out -e task3.err
#SBATCH --gres=gpu:1

if [ "$1" = "compile" ]; then
  	module load cuda
  	nvcc task3.cu vscale.cu -Xcompiler -O3 -Xcompiler -Wall -Xptxas -O3 -o task3
elif [ "$1" = "run" ]; then
  	for i in {10..29}; do
    	n=$((2 ** i))
    	echo "Running with n = 2^$i = $n"
    	./task3 $n
    	echo
  	done
elif [ "$1" = "clean" ]; then
  	rm -f task3 task3.out task3.err
else
  	echo "./$0 [compile | run | clean]"
fi
