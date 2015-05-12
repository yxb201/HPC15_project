#include "p3dfft.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{	
	int nproc, proc_id; 
	int dims[2], memsize[3];
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id)

	// 2x2 processor for now
	dims[0]=2; dims[1]=2;

	printf("pi = %f\n",M_PI);
	
	MPI_Finalize();


}
