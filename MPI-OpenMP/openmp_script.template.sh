#!/bin/bash
#
# One node is set for execution
#SBATCH -N 1 
# Only one process is needed for OpenMP intra-node parallelism
#SBATCH --ntasks=1 
# Select the maximum number of OpenMP threads that you will use
#SBATCH --cpus-per-task=2
#SBATCH --partition=rome
#
# Five minutes are set for execution
#SBATCH -t 00:05:00
#

# select the number of OpenMP threads for the given execution
export OMP_NUM_THREADS=2
# running MPI job with 2 cores
./a.out
