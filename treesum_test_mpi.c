
/*
 * File:     treesum_test_mpi.c
 *
 * Purpose:  test tree-structured communication to find the global sum
 *           of a random collection of ints.  This version doesn't
 *           require that no of processes  be a power of 2.
 *
 * Compile:  mpicc -Wall -o treesum_test_mpi  treesum_test_mpi.c treesum_mpi.c
 * minunit.c Run:      ibrun  -v ./treesum_test_mpi
 *
 * IPP:      Programming Assignment Text Book Problem 3.3 Page 148
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "minunit.h"

int my_rank, no_proc;
MPI_Comm comm;
char msgbuf[100];

int global_sum(int my_int, int my_rank, int no_proc, MPI_Comm comm);

const int nums[16] = { 4,7,3,9,19,35,-3,98, 9, -7, 8, 1, 13,15,14,16};

/*-------------------------------------------------------------------
 * Test global_sum()
 * If successful, return NULL
 */
char *treesum_test1() {
    if (no_proc >= 16){ 
        printf("too many processors\n");
        exit(16);
    }

    double startwtime = 0, endwtime = 0;
    double sum, expected = 0;

    for (int i = 0; i < no_proc; i++){ 
        expected += nums[i]; 
    }

    sum = global_sum(nums[my_rank], my_rank, no_proc, comm);
    if (my_rank == 0) startwtime = MPI_Wtime();
    if(my_rank == 0){ 
        endwtime =  MPI_Wtime();
        printf("%s: Wall clock time = %f at Proc 0 of %d processes\n", "Test 1",
           endwtime - startwtime, no_proc);
    printf("%s: sum  = %f\n", "Test 1", sum);
    printf("%s: expected  = %f\n", "Test 1", expected);
    sprintf(msgbuf, "%s: Failed due to a large error", "Test 1");
    mu_assert(msgbuf, sum == expected);
    }

    return NULL;
}

/*-------------------------------------------------------------------
 * Run all tests.  Ignore returned messages.
 */
void run_all_tests(void) { mu_run_test(treesum_test1); }

/*-------------------------------------------------------------------
 * The main entrance to run all tests.
 * Only Proc 0 prints the test summary
 */
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &no_proc);
  MPI_Comm_rank(comm, &my_rank);

  run_all_tests();

  if (my_rank == 0) {
    mu_print_test_summary("Summary:");
  }
  MPI_Finalize();
  return 0;
}
