#include "p3dfft.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{	
	int nproc, proc_id, conf; 
	int dims[2], memsize[3];
	double L = 2.0 * M_PI; 
	int nx,ny,nz;
	int istart[3], isize[3], iend[3];
	int fstart[3], fsize[3], fend[3];


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

	// 2x2 processor for now
	dims[0]=2; dims[1]=2;
	
	// assume for now 
	nx = ny = nz = 4;
	
	// initialize P3DFFT
	Cp3dfft_setup(dims,nx,ny,nz,MPI_Comm_c2f(MPI_COMM_WORLD), nx,ny,nz, 0, memsize);
	
	conf = 1;
	Cp3dfft_get_dims(istart, iend, isize, conf);

if (proc_id == 3){
	printf("istart[0] = %d, istart[1] = %d, istart[2] = %d\n", istart[0],istart[1], istart[2]);
	printf("iend[0] = %d, iend[1] = %d, iend[2] = %d\n", iend[0],iend[1], iend[2]);
	printf("isize[0] = %d, isize[1] = %d, isize[2] = %d\n", isize[0],isize[1], isize[2]);
}

	
	Cp3dfft_clean();


	MPI_Finalize();


}
