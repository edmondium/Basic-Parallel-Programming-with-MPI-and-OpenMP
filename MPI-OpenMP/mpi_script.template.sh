#!/bin/bash
#
# One node is set for execution
#SBATCH -N 1 
# Select the number of MPI processes per node
# #SBATCH --ntasks=XXX
#SBATCH --partition=rome
#
# Five minutes are set for execution
#SBATCH -t 00:05:00
#
# running MPI job with 2 cores
srun -n 2 ./a.out
