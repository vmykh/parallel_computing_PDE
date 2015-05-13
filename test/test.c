#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char* argv[])
{
  printf("before init\n");


  MPI_Init (&argc, &argv);
	
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int PROCS_AMOUNT;
  MPI_Comm_size(MPI_COMM_WORLD, &PROCS_AMOUNT);

  printf("after init\n");

  MPI_Finalize();

  printf("after finalize\n");
  
  return 0;
}
