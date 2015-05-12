#include "p3dfft.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{	
	int nproc, proc_id, proc_i, proc_j, conf; 
	int dims[2], memsize[3];
	double L, h; 
	int nx,ny,nz;
	int istart[3], isize[3], iend[3];
	int fstart[3], fsize[3], fend[3];
	int *spread_data;	
	int Msp = 2; 
	int west, east, north, south, ne, nw, se, sw, P1, P2;
	int i,j,k;
	L = 2.0 * M_PI;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

	// 2x2 processor for now
	dims[0]=3; dims[1]=4;
	P1 = dims[0]; P2 = dims[1];	

	int *proc_grid = (int *) malloc( sizeof(int)*(P1+2)*(P2+2) );
	// processor 0: generate a grid of processor IDs (with neighbours padded)
	if (proc_id == 0){

		for(i=1; i<=P2; ++i)
			for(j=1; j<=P1; ++j)
				proc_grid[i*(P1+2)+j] = (i-1)*P1 + j-1;

		// west and east boundaries 
		for(i=1; i<=P2; ++i){
			proc_grid[i*(P1+2)] = proc_grid[i*(P1+2)+P1]; 
			proc_grid[i*(P1+2)+P1+1] = proc_grid[i*(P1+2)+1];	  
		}

		// north and south boundaries
		for(j=1; j<=P1; ++j ){
			proc_grid[j] = proc_grid[ P2*(P1+2) + j];
			proc_grid[(P2+1)*(P1+2)+j] = proc_grid[(P1+2)+j];
		}

		// corners: sw, ne, se, nw
		proc_grid[0] = proc_grid[(P1+2)*P2+P1];
		proc_grid[(P2+2)*(P1+2)-1] = proc_grid[P1+2 + 1];
		proc_grid[P1+1] = proc_grid[(P1+2)*P2+1];
		proc_grid[(P1+2)*(P2+1)] = proc_grid[(P1+2)+P1];

		for(i=0; i<P2+2; ++i){
			for (j=0; j<P1+2; ++j){
				printf("%2d ", proc_grid[i*(P1+2)+j]);}
			printf("\n");	
		}
	
	}
	
	MPI_Bcast(proc_grid, (P1+2)*(P2+2), MPI_INT, 0, MPI_COMM_WORLD );

	/*
	for(i=0; i<P2+2; ++i){
		for (j=0; j<P1+2; ++j){
			printf("%2d ", proc_grid[i*(P1+2)+j]);}
		printf("\n");	
	}
	*/

	// find the indices of proc_id  in the proc_grid
	proc_i = proc_id / P1 + 1;
	proc_j = proc_id % P1 + 1;
	
	// set neighbours
	west  = proc_grid[  proc_i   *(P1+2) + proc_j-1 ];
	east  = proc_grid[  proc_i   *(P1+2) + proc_j+1 ];
	north = proc_grid[ (proc_i+1)*(P1+2) + proc_j   ];
	south = proc_grid[ (proc_i-1)*(P1+2) + proc_j   ];
	nw    = proc_grid[ (proc_i+1)*(P1+2) + proc_j-1 ];
	ne    = proc_grid[ (proc_i+1)*(P1+2) + proc_j+1 ];
	sw    = proc_grid[ (proc_i-1)*(P1+2) + proc_j-1 ];
	se    = proc_grid[ (proc_i-1)*(P1+2) + proc_j+1 ];

//	printf("w  = %2d, e  = %2d, n  = %2d, s  = %2d\n", west, east, north, south);
//	printf("nw = %2d, ne = %2d, sw = %2d, se = %2d\n", nw, ne, sw, se );

	

	free(proc_grid);


	// assume for now 
	nx = ny = nz = 12;

	h = L / nx;
	
	// initialize P3DFFT
	Cp3dfft_setup(dims,nx,ny,nz,MPI_Comm_c2f(MPI_COMM_WORLD), nx,ny,nz, 0, memsize);
	
	conf = 1;
	Cp3dfft_get_dims(istart, iend, isize, conf);


	// allocate spread data
	spread_data = (int *) malloc(sizeof(int) * isize[0]*isize[1]*isize[2]);

if (proc_id == 3){

	
	printf("istart[0] = %d, istart[1] = %d, istart[2] = %d\n", istart[0],istart[1], istart[2]);
	printf("iend[0] = %d, iend[1] = %d, iend[2] = %d\n", iend[0],iend[1], iend[2]);
	printf("isize[0] = %d, isize[1] = %d, isize[2] = %d\n", isize[0],isize[1], isize[2]);
	
	// generate a source 
	double xj[3];
	xj[0]=6.5*h; xj[1]=6.5*h; xj[2]=6.5*h;
	int xi[3];
	xi[0] = (int) (xj[0]/h); 
	xi[1] = (int) (xj[1]/h);
	xi[2] = (int) (xj[2]/h);
	
	printf("%.16f \n", xj[0]/h);
	printf("%d %d %d \n", xi[0], xi[1], xi[2] );
}
	Cp3dfft_clean();
	
	free(spread_data);

	MPI_Finalize();


}
