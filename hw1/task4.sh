#!/usr/bin/env bash

#SBATCH -p wacc

#SBATCH -c 2
#SBATCH -J FirstSlurm
#SBATCH -o FirstSlurm.out -e FirstSlurm.err
hostname
