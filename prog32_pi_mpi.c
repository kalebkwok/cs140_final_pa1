
/*
 * File:     prog32_pi_mpi.c
 *
 * Purpose:  Estimate pi using a monte carlo method
 *
 * IPP book:      Programming Assignment 3.2
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/*-------------------------------------------------------------------
 * Estimate Pi in parallel using a monte carlo method
 * Only process 0 returns the correct pi estimation
 * In args:
 * 	number_of_tosses: No of tosses.
 *      my_rank: my process rank
 *      no_proc: number of processes
 *      comm: MPI communication group
 * Return: estimated pi
 */

double parallel_pi(long long int number_of_tosses, int my_rank, int no_proc,
                   MPI_Comm comm) {
    long long int i;
    double x, y,distance_squared;
    double pi_estimate = 0.0;
    long long int number_in_circle = 0;
    
    for (i = 0; i < number_of_tosses / (long long int) no_proc; i++) {
      x = 2 * random() / ((double)RAND_MAX) - 1.0;
      y = 2 * random() / ((double)RAND_MAX) - 1.0;
      distance_squared = x * x + y * y;
      if (distance_squared <= 1) {
        number_in_circle++;
      }
    }
    
    long long int* total_hit = malloc(sizeof(long long int));
    *total_hit = 0;
    
    MPI_Reduce(&number_in_circle, total_hit, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, comm);
    
    if(my_rank == 0){
      pi_estimate = 4*((double)(*total_hit))/((double)number_of_tosses);
    }
    
    free(total_hit);

    return pi_estimate;
}
