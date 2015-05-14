#include "p3dfft.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define l(i,j,k,nx,ny) ((i)+(j)*(nx)+(k)*(nx)*(ny)) 

int mod(int a, int b){
	
	int r = a % b; 
	return r < 0 ? r+b : r;
	
}

int main(int argc, char *argv[])
{	
	int nproc, proc_id, proc_i, proc_j, conf; 
	int dims[2], memsize[3];
	double L, h; 
	int nx,ny,nz;
	int lnx, lny, lnz;
	int istart[3], isize[3], iend[3];
	int fstart[3], fsize[3], fend[3];
	double *spread_rect, *local_rect;	
	int Msp = 2; 
	int mx, my, mz, lmx, lmy, lmz, smx, smy, smz;
	int NORTH, NE, EAST, SE, SOUTH, SW, WEST, NW, P1, P2;
	int i,j,k;
 	double *N_Recv, *NE_Recv, *E_Recv, *SE_Recv, *S_Recv, *SW_Recv, *W_Recv, *NW_Recv;
	double *N_Send, *NE_Send, *E_Send, *SE_Send, *S_Send, *SW_Send, *W_Send, *NW_Send;	

	L = 2.0 * M_PI;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Status status;

	// 2x2 processor for now
	dims[0]=3; dims[1]=4;
	P1 = dims[0]; P2 = dims[1];	

	// set 8 neighbours
	k = proc_id / P1;
	j = proc_id % P1;
	NORTH = mod(k+1,P2)*P1 + mod(j  ,P1);
	NE    = mod(k+1,P2)*P1 + mod(j+1,P1);
	EAST  = mod(k  ,P2)*P1 + mod(j+1,P1);
	SE    = mod(k-1,P2)*P1 + mod(j+1,P1);
	SOUTH = mod(k-1,P2)*P1 + mod(j  ,P1);
	SW    = mod(k-1,P2)*P1 + mod(j-1,P1);
	WEST  = mod(k  ,P2)*P1 + mod(j-1,P1);
	NW	  = mod(k+1,P2)*P1 + mod(j-1,P1);
/* 
if (proc_id == 11){
	printf("N=%2d, NE=%2d, E=%2d, SE=%2d, S =%2d, SW=%2d, W=%2d, NW=%2d\n", NORTH,NE,EAST,SE,SOUTH,SW,WEST,NW);
}
*/

	// assume for now 
	nx = ny = nz = 12;

	h = L / nx;
	
	// initialize P3DFFT
	Cp3dfft_setup(dims,nx,ny,nz,MPI_Comm_c2f(MPI_COMM_WORLD), nx,ny,nz, 0, memsize);
	
	conf = 1;
	Cp3dfft_get_dims(istart, iend, isize, conf);

	lnx = isize[0];
	lny = isize[1];
	lnz = isize[2];

	// rectangle for spreading, dimension: nx x (ny_local + 2Msp) x (nz_local+2Msp)
	spread_rect = (double *) malloc( sizeof(double) * lnx*(lny+2*Msp)*(lnz+2*Msp) );
	for (i=0; i<lnx*(lny+2*Msp)*(lnz+2*Msp); ++i) spread_rect[i] = 0.;

	// rectangle for local data, dimension: nx x ny_local x nz_local
	local_rect  = (double *) malloc( sizeof(double) * lnx*lny*lnz );	
	for (i=0; i<lnx*lny*lnz; ++i) local_rect[i] = 0.;


	// allocate buffer size
	N_Recv = (double *) malloc( sizeof(double) * lnx*lny*Msp );
	N_Send = (double *) malloc( sizeof(double) * lnx*lny*Msp );
	S_Recv = (double *) malloc( sizeof(double) * lnx*lny*Msp );
	S_Send = (double *) malloc( sizeof(double) * lnx*lny*Msp );
	W_Recv = (double *) malloc( sizeof(double) * lnx*Msp*lnz );
	W_Send = (double *) malloc( sizeof(double) * lnx*Msp*lnz );
	E_Recv = (double *) malloc( sizeof(double) * lnx*Msp*lnz );
	E_Send = (double *) malloc( sizeof(double) * lnx*Msp*lnz );

	NW_Recv = (double *) malloc( sizeof(double) * lnx*Msp*Msp );
	NW_Send = (double *) malloc( sizeof(double) * lnx*Msp*Msp );
	SW_Recv = (double *) malloc( sizeof(double) * lnx*Msp*Msp );
	SW_Send = (double *) malloc( sizeof(double) * lnx*Msp*Msp );
	NE_Recv = (double *) malloc( sizeof(double) * lnx*Msp*Msp );
	NE_Send = (double *) malloc( sizeof(double) * lnx*Msp*Msp );
	SE_Recv = (double *) malloc( sizeof(double) * lnx*Msp*Msp );
	SE_Send = (double *) malloc( sizeof(double) * lnx*Msp*Msp );

	
if (proc_id == 7){

	
	printf("istart[0] = %d, istart[1] = %d, istart[2] = %d\n", istart[0],istart[1], istart[2]);
	printf("iend[0] = %d, iend[1] = %d, iend[2] = %d\n", iend[0],iend[1], iend[2]);
	printf("isize[0] = %d, isize[1] = %d, isize[2] = %d\n", isize[0],isize[1], isize[2]);
	
	// generate a source 
	double xj[3];
	xj[0]=6.5*h; xj[1]=6.5*h; xj[2]=6.5*h;
	
	mx = (int) (xj[0]/h); 
	my = (int) (xj[1]/h);
	mz = (int) (xj[2]/h);
	
	printf("center: %d %d %d \n", mx, my, mz);

	// mx,my,mz in the spread_rect
	smx= my - (istart[0]-1);
	smy= my - (istart[1]-1) + Msp;
	smz= mz - (istart[2]-1) + Msp;

	printf("lnx = %d, lny = %d, lnz = %d \n", lnx, lny, lnz);
	printf("smx = %d, smy = %d, smz = %d \n", smx, smy, smz);

	for (k=-Msp+1; k<=Msp; ++k){
	for (j=-Msp+1; j<=Msp; ++j){
	for (i=-Msp+1; i<=Msp; ++i){ 

//		printf("(i,j,k) = (%d, %d, %d) \n", i+smx, j+smy, k+smz);
		spread_rect[l(mod(i+smx,nx),j+smy,k+smz,nx,lny+2*Msp)] += 1.;
	}
	}
	}

	for (k=0; k < Msp; ++k){
		for (j=0; j<Msp; ++j){
			for (i=0; i<lnx; ++i){
				SE_Send[l(i,j,k,nx,Msp)] = spread_rect[l(0+i,Msp+lny+j,0+k,nx, lny+2*Msp)];
				printf("%1.1f ", SE_Send[l(i,j,k,nx,Msp)]);
			}
			printf("\n");
		}
		printf("\n\n\n");
	}

	MPI_Send(SE_Send, lnx*Msp*Msp, MPI_DOUBLE, SE, 99, MPI_COMM_WORLD);

/*
	
	for (k=1; k<=4; ++k){
	for (j=3; j<=6; ++j){
	for (i=5; i<=8; ++i){

//		spread_rect[l(i,j,k,nx,lny+2*Msp)] = 1;		
//		printf("(i,j,k) = (%d, %d, %d) \n", i, j, k);
		printf("%1.2f ", spread_rect[l(i,j,k,nx,lny+2*Msp)]);
//		printf("%1.2f ", spread_rect[i + nx*j + nx*(lny+2*Msp)*k]);
	}
		printf("\n");
	}
		printf("\n \n \n");
	}

//	double tmp = 0;
	for (k=0; k<lnz+2*Msp; ++k){
	for (j=0; j<lny+2*Msp; ++j){
	for (i=0; i<lnx ; ++i){
//		tmp = tmp + spread_rect[l(i,j,k,nx,lny+2*Msp)];
		printf("%1.1f ", spread_rect[l(i,j,k,nx,lny+2*Msp)]);
	}
		printf("\n");
	}
		printf("\n \n \n ");
	}
//	printf("final sum = %f\n", tmp);	

*/
}


if (proc_id == 5){
	
	printf("NW=%d\n",NW);
	MPI_Recv(NW_Recv, lnx*Msp*Msp, MPI_DOUBLE, NW, 99, MPI_COMM_WORLD, &status);
	for(i=0; i<lnx*Msp*Msp; ++i)
		printf("%1.1f ", NW_Recv[i]);	
	
	printf("received\n");

}


	Cp3dfft_clean();
	
	free(spread_rect);
	free(local_rect);
	free(N_Recv); free(N_Send); free(S_Recv); free(S_Send); free(W_Recv); free(W_Send);  free(E_Recv); free(E_Send); 
	free(NW_Recv); free(NW_Send); free(SE_Recv); free(SE_Send); free(SW_Recv); free(SW_Send);  free(NE_Recv); free(NE_Send); 


	MPI_Finalize();

	return 0;
}
