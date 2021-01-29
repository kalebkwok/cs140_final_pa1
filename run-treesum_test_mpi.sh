#!/bin/bash
#SBATCH --job-name="treesum"
#SBATCH --output="job_treesum4.%j.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --export=ALL
#SBATCH -t 00:02:00

#This job runs with 1 node, 4 cores per node for a total of 4 cores.
#Total number of MPI processes= 1*4=4
#ibrun in verbose mode will give binding detail

ibrun -v ./treesum_test_mpi
