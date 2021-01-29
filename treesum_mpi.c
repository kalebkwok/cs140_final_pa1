/*
 * File:     treesum_mpi.c
 *
 * Purpose:  Use tree-structured communication to find the global sum
 *           of a random collection of ints.  This version doesn't
 *           require that comm_sz be a power of 2.
 *
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/*-------------------------------------------------------------------
 * Function:
 *  global_sum
 *
 * Purpose:
 *  Implement a global sum using tree-structured communication
 *
 * Notes:
 *  1.  The return value for global sum is only valid on process 0
 */
int global_sum(int my_int /* in */, int my_rank /* in */, int comm_sz /* in */,
               MPI_Comm comm /* in */) {
    
    int my_sum = my_int;

    for(int i = 1; i<comm_sz; i*=2){
        if(my_rank % (i * 2) == 0){
            if(my_rank + i < comm_sz){
                int buff;
                MPI_Recv(&buff, 1, MPI_INT, my_rank + i, my_rank, comm, MPI_STATUS_IGNORE);
                my_sum+=buff;
            }
        }
        else{
            MPI_Send(&my_sum, 1, MPI_INT, my_rank - i, my_rank - i, comm);
            return my_sum;
        }
    }
   return my_sum;
} /* Global_sum */
