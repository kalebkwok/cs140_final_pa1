#!/bin/bash
#Next line shows the job name you can find when querying the job status
#SBATCH --job-name="itmv"
#Next line is the output file name of the execution log
#SBATCH --output="job_itmv_2nodes.%j.out"
#Next line shows where to ask for machine nodes
#SBATCH --partition=compute
#Next line asks for 1 node and  2 cores per node for a total of 2 cores.
#Total number of MPI processes= 1*2=2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --export=ALL
#Next line limits the job execution time at most 3 minute.
#SBATCH -t 00:03:00

#ibrun in verbose mode will give binding detail

ibrun -v ./itmv_mult_test_mpi
