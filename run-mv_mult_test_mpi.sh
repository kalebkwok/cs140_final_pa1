#!/bin/bash
#Next line shows the job name you can find when querying the job status
#SBATCH --job-name="mv_mpi"
#Next line is the output file name of the execution log
#SBATCH --output="job_mv_4nodes.%j.out"
#Next line shows where to ask for machine nodes
#SBATCH --partition=compute
#This job runs with 1 node, 4 cores per node for a total of 4 cores.
#Total number of MPI processes= 1*4=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --export=ALL
#Next line limits the job execution time at most 3 minutes.
#SBATCH -t 00:03:00
#ibrun in verbose mode will give binding detail

ibrun -v ./mv_mult_test_mpi
