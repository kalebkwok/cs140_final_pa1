#!/bin/bash
#Next line shows the job name you can find when querying the job status
#SBATCH --job-name="prog32pi_mpi"
#Next line is the output file name of the execution log
#SBATCH --output="job_prog32_pi2.%j.out"
#Next line shows where to ask for machine nodes
#SBATCH --partition=compute
#Next line asks for 1 node and  4 cores per node for a total of 4 cores.
#Total number of MPI processes= 1*4=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --export=ALL
#Next line limits the job execution time at most 3 minute.
#SBATCH -t 00:02:00

#ibrun in verbose mode will give binding detail

ibrun -v ./prog32_pi_test_mpi
