
/*
 * File:     prog32_pi_test_mpi.c
 *
 * Purpose:  Estimate pi using a monte carlo method
 *
 * IPP:      Programming Assignment 3.2 (Text Book Page 148)
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "minunit.h"

#include <math.h>

double PI25DT = 3.141592653589793238462643;
int my_rank, no_proc;
MPI_Comm comm;
char msgbuf[100];

double parallel_pi(long long int number_of_tosses, int my_rank, int no_proc,
                   MPI_Comm comm);

/*-------------------------------------------------------------------
 * Test with a number  of tosses
 * Only Process 0 conducts a result compiarsion. If failed, return a message
 * string If successful, return NULL
 */
char *pi_test(char *testmsg, long long int number_of_tosses, double threshold) {
  double startwtime = 0, endwtime = 0;
  double error, pi_estimate;

  if (my_rank == 0) startwtime = MPI_Wtime();

  srandom(my_rank + 1);
  pi_estimate = parallel_pi(number_of_tosses, my_rank, no_proc, comm);
  if (my_rank == 0) {
    endwtime = MPI_Wtime();
    printf("%s: Wall clock time = %f at Proc 0 of %d processes\n", testmsg,
           endwtime - startwtime, no_proc);
    printf("%s: pi estimate = %f\n", testmsg, pi_estimate);
    error = fabs(pi_estimate - PI25DT);
    printf("%s: Error = %f\n", testmsg, error);
    sprintf(msgbuf, "%s: Failed due to a large error", testmsg);
    mu_assert(msgbuf, error < threshold);
  }
  return NULL;
}

char *pi_test1() { return pi_test("Test 1", 1000000, 0.01); }

char *pi_test2() { return pi_test("Test 2", 100000000, 0.001); }

char *pi_test3() { return pi_test("Test 3", 1000000000, 0.0001); }

/*-------------------------------------------------------------------
 * Run all tests.  Ignore returned messages.
 */
void run_all_tests(void) {
  mu_run_test(pi_test1);
  mu_run_test(pi_test2);
  mu_run_test(pi_test3);
}

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
