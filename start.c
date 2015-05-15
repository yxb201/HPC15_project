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


void setSbuffer(double *buffer, double *rect, int x[3], int dimbuffer[3], int dimrect[3] ){
	
	int i,j,k;
	
	for (k=0; k<dimbuffer[2]; ++k){
		for (j=0; j<dimbuffer[1]; ++j){
			for (i=0; i<dimbuffer[0]; ++i){
				buffer[l(i,j,k,dimbuffer[0], dimbuffer[1])] = rect[l(x[0]+i,x[1]+j,x[2]+k, dimrect[0], dimrect[1])];	
			}
		}
	}

}


void getRbuffer(double *buffer, double *rect, int x[3], int dimbuffer[3], int dimrect[3]){
	int i,j,k;

	for (k=0; k<dimbuffer[2]; ++k){
		for (j=0; j<dimbuffer[1]; ++j){
			for (i=0; i<dimbuffer[0]; ++i){
				rect[l(x[0]+i,x[1]+j,x[2]+k,dimrect[0],dimrect[1])] += buffer[l(i,j,k,dimbuffer[0], dimbuffer[1])];
			}
		}
	}
}


int main(int argc, char *argv[])
{	
	int nproc, proc_id, proc_i, proc_j, conf; 
	int dims[2], memsize[3];
	double L, h, tau, diffx, diffy, diffz, E1, E2x, E2y, E2z; 
	int M, Mr, R;
	int lnx, lny, lnz;
	int istart[3], isize[3], iend[3];
	int fstart[3], fsize[3], fend[3];
	double *spread_rect, *local_rect;
	double *E2xl, *E2yl, *E2zl;	
	int Msp = 2; 
	int mx, my, mz, lmx, lmy, lmz, smx, smy, smz;
	int NORTH, NE, EAST, SE, SOUTH, SW, WEST, NW, P1, P2;
	int i,j,k,s, l1, l2, l3;
 	double *N_Recv, *NE_Recv, *E_Recv, *SE_Recv, *S_Recv, *SW_Recv, *W_Recv, *NW_Recv;
	double *N_Send, *NE_Send, *E_Send, *SE_Send, *S_Send, *SW_Send, *W_Send, *NW_Send;	
	int dimSbuffer[3], dimRbuffer[3], dimSpreadRect[3];
	int idx[3]; 
	L = 2.0 * M_PI;
	double *xj, *yj, *zj;
	int n_src, N_src;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Status status;

	// 2x2 processor for now
	dims[0]=3; dims[1]=4;
	P1 = dims[0]; P2 = dims[1];	

	// number of sources in each processor
	n_src = 1;

	// number of sources in total
	N_src = n_src * P1 * P2;


	xj = (double *) malloc( sizeof(double) * n_src );
	yj = (double *) malloc( sizeof(double) * n_src );
	zj = (double *) malloc( sizeof(double) * n_src );

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

	// assume for now 
	M   = 12; 
	R   = 2;
	Mr  = R*M; 
	tau = (1.*Msp) / (M*M); 

	h = L / M;
	
	E2xl = (double *) malloc( sizeof(double) * 2*Msp );
	E2yl = (double *) malloc( sizeof(double) * 2*Msp );
	E2zl = (double *) malloc( sizeof(double) * 2*Msp );

	// initialize P3DFFT
	Cp3dfft_setup(dims,M,M,M,MPI_Comm_c2f(MPI_COMM_WORLD), M,M,M, 0, memsize);
	
	conf = 1;
	Cp3dfft_get_dims(istart, iend, isize, conf);

	// dimension of local rectangle
	lnx = isize[0];
	lny = isize[1];
	lnz = isize[2];
	
	// one source for now
	xj[0]= (istart[0]+iend[0]*)h/2; 
	yj[0]= (istart[1]+iend[1])*h/2; 
	xj[0]= (istart[2]+iend[2])*h/2;
	

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



//if (proc_id == 9){
/*
	printf("istart[0] = %d, istart[1] = %d, istart[2] = %d\n", istart[0],istart[1], istart[2]);
	printf("iend[0] = %d, iend[1] = %d, iend[2] = %d\n", iend[0],iend[1], iend[2]);
	printf("isize[0] = %d, isize[1] = %d, isize[2] = %d\n", isize[0],isize[1], isize[2]);
*/

	// for each source
	double mxh  = mx*h;
	double myh  = my*h;
	double mzh  = mz*h;
	double piMtau = M_PI / (Mr * tau);
	for(s=0; s < n_src ; ++s){
		
		// find the closest grid point (in the whole domain)
		mx = (int) ( xj[s]/h );
		my = (int) ( yj[s]/h );
		mz = (int) ( zj[s]/h );
		
		// closest grid point (in spreading rect with halo cells )
		smx= mx - (istart[0]-1);
		smy= my - (istart[1]-1) + Msp;
		smz= mz - (istart[2]-1) + Msp;

		diffx = xj[s] - mxh;
		diffy = yj[s] - myh;
		diffz = zj[s] - mzh;
		E1 = exp( -(diffx*diffx+diffy*diffy+diffz*diffz)/(4*tau) );

		E2x = exp( piMtau * diffx  );
		E2y = exp( piMtau * diffy  );
		E2z = exp( piMtau * diffz  );

		for (l1 = -Msp+1; l1<=Msp; ++l1 ){
			//////////
		}


	}

/* debug message 
if (proc_id == 11){	
	printf("center: %d %d %d \n", mx, my, mz);
}
*/
	// mx,my,mz in the spread_rect
//	printf("lnx = %d, lny = %d, lnz = %d \n", lnx, lny, lnz);
//	printf("smx = %d, smy = %d, smz = %d \n", smx, smy, smz);

	dimSpreadRect[0]=lnx; dimSpreadRect[1]=lny+2*Msp; dimSpreadRect[2] = lnz+2*Msp;

	// build the spreading rectangle
	for (k=-Msp+1; k<=Msp; ++k){
		for (j=-Msp+1; j<=Msp; ++j){
			for (i=-Msp+1; i<=Msp; ++i){ 
				spread_rect[l(mod(i+smx,nx),j+smy,k+smz,lnx,lny+2*Msp)] += 1.;
			}
		}
	}

//}


	// copy spreading rectangle to local rectangle
	idx[0] = 0; idx[1] = Msp; idx[2] = Msp;
	setSbuffer(local_rect, spread_rect, idx, isize, dimSpreadRect );


	// set North Send buffer
	idx[0] = 0; idx[1]=Msp; idx[2]=lnz+Msp;
	dimSbuffer[0] = lnx; dimSbuffer[1]=lny; dimSbuffer[2]=Msp; 	
	setSbuffer(N_Send, spread_rect, idx, dimSbuffer, dimSpreadRect);

	// set South Send buffer
	idx[0] = 0; idx[1]=Msp; idx[2]=0; 	
	setSbuffer(S_Send, spread_rect, idx, dimSbuffer, dimSpreadRect);

	// set West Send buffer
	idx[0] = 0; idx[1]=0; idx[2]=Msp;
	dimSbuffer[0] = lnx; dimSbuffer[1]=Msp; dimSbuffer[2]=lnz; 	
	setSbuffer(W_Send, spread_rect, idx, dimSbuffer, dimSpreadRect);
	
	// set East Send buffer
	idx[0] = 0; idx[1]=lny+Msp; idx[2]=Msp; 	
	setSbuffer(E_Send, spread_rect, idx, dimSbuffer, dimSpreadRect);

	// set NE Send buffer
	idx[0] = 0; idx[1]=lny+Msp; idx[2]=lnz+Msp;
	dimSbuffer[0] = lnx; dimSbuffer[1]=Msp; dimSbuffer[2]=Msp; 	
	setSbuffer(NE_Send, spread_rect, idx, dimSbuffer, dimSpreadRect);

	// set SE send buffer
	idx[0] = 0; idx[1]=lny+Msp; idx[2]=0; 	
	setSbuffer(SE_Send, spread_rect, idx, dimSbuffer, dimSpreadRect);

	// set SW send buffer
	idx[0] = 0; idx[1]=0; idx[2]=0; 	
	setSbuffer(SW_Send, spread_rect, idx, dimSbuffer, dimSpreadRect);

	// set NW send buffer
	idx[0] = 0; idx[1]=0; idx[2]=lnz+Msp; 	
	setSbuffer(NW_Send, spread_rect, idx, dimSbuffer, dimSpreadRect);


/*
	
	printf("center: %d %d %d \n", mx, my, mz);

	for(k=0; k<lnz; ++k){
		for(j=0; j<lny; ++j){
			for(i=0; i<lnx; ++i){
				printf("%1.1f ", local_rect[l(i,j,k,lnx,lny)]);
			}
			printf("\n");
		}
		printf("\n\n");
	}
*/

/*
	for (k=0; k < Msp; ++k){
		for (j=0; j<Msp; ++j){
			for (i=0; i<lnx; ++i){
//				SE_Send[l(i,j,k,nx,Msp)] = spread_rect[l(0+i,Msp+lny+j,0+k,nx, lny+2*Msp)];
				printf("%1.1f ", SE_Send[l(i,j,k,nx,Msp)]);
			}
			printf("\n");
		}
		printf("\n\n\n");
	}
*/


	// MPI send to neighbours
	MPI_Send( N_Send, lnx*lny*Msp, MPI_DOUBLE, NORTH, 1, MPI_COMM_WORLD);
	MPI_Send( S_Send, lnx*lny*Msp, MPI_DOUBLE, SOUTH, 2, MPI_COMM_WORLD);
	MPI_Send( W_Send, lnx*Msp*lnz, MPI_DOUBLE,  WEST, 3, MPI_COMM_WORLD);
	MPI_Send( E_Send, lnx*Msp*lnz, MPI_DOUBLE,  EAST, 4, MPI_COMM_WORLD);
	MPI_Send(NE_Send, lnx*Msp*Msp, MPI_DOUBLE,    NE, 5, MPI_COMM_WORLD);
	MPI_Send(SE_Send, lnx*Msp*Msp, MPI_DOUBLE,    SE, 6, MPI_COMM_WORLD);
	MPI_Send(SW_Send, lnx*Msp*Msp, MPI_DOUBLE,    SW, 7, MPI_COMM_WORLD);
	MPI_Send(NW_Send, lnx*Msp*Msp, MPI_DOUBLE,    NW, 8, MPI_COMM_WORLD);


	// MPI receive from neighbours
	MPI_Recv( N_Recv, lnx*lny*Msp, MPI_DOUBLE, NORTH, 2, MPI_COMM_WORLD, &status);
	MPI_Recv( S_Recv, lnx*lny*Msp, MPI_DOUBLE, SOUTH, 1, MPI_COMM_WORLD, &status);
	MPI_Recv( W_Recv, lnx*Msp*lnz, MPI_DOUBLE,  WEST, 4, MPI_COMM_WORLD, &status);
	MPI_Recv( E_Recv, lnx*Msp*lnz, MPI_DOUBLE,  EAST, 3, MPI_COMM_WORLD, &status);
	MPI_Recv(NE_Recv, lnx*Msp*Msp, MPI_DOUBLE,    NE, 7, MPI_COMM_WORLD, &status);
	MPI_Recv(SE_Recv, lnx*Msp*Msp, MPI_DOUBLE,    SE, 8, MPI_COMM_WORLD, &status);
	MPI_Recv(SW_Recv, lnx*Msp*Msp, MPI_DOUBLE,    SW, 5, MPI_COMM_WORLD, &status);
	MPI_Recv(NW_Recv, lnx*Msp*Msp, MPI_DOUBLE,    NW, 6, MPI_COMM_WORLD, &status);


	// copy receive buffer to local rectangle
	// add contribution from N buffer
	idx[0]=0; idx[1]=0; idx[2]=lnz-Msp;
	dimRbuffer[0]=lnx; dimRbuffer[1]=lny; dimRbuffer[2]=Msp; 
	getRbuffer( N_Recv, local_rect, idx, dimRbuffer, isize );

	// add contribution from S buffer
	idx[0]=0; idx[1]=0; idx[2]=0;
	getRbuffer( S_Recv, local_rect, idx, dimRbuffer, isize );

	// add contribution from W buffer
	idx[0]=0; idx[1]=0; idx[2]=0;
	dimRbuffer[0]=lnx; dimRbuffer[1]=Msp; dimRbuffer[2]=lnz;
	getRbuffer( W_Recv, local_rect, idx, dimRbuffer, isize );
	
	// add contribution from E buffer
	idx[0]=0; idx[1]=lny-Msp; idx[2]=0;
	getRbuffer( E_Recv, local_rect, idx, dimRbuffer, isize );

	// add contribution form NW buffer 
	idx[0]=0; idx[1]=0; idx[2]=lnz-Msp;
	dimRbuffer[0]=lnx; dimRbuffer[1]=Msp; dimRbuffer[2]=Msp; 
	getRbuffer( NW_Recv, local_rect, idx, dimRbuffer, isize );

	// add contribution from SW buffer 
	idx[0]=0; idx[1]=0; idx[2]=0;
	getRbuffer( SW_Recv, local_rect, idx, dimRbuffer, isize );

	// add contribution from SE buffer
	idx[0]=0; idx[1]=lny-Msp; idx[2]=0;
	getRbuffer( SE_Recv, local_rect, idx, dimRbuffer, isize );

	// add contribution from NE buffer
	idx[0]=0; idx[1]=lny-Msp; idx[2]=lnz-Msp;
	getRbuffer( NE_Recv, local_rect, idx, dimRbuffer, isize );



if (proc_id == 2){ 

	for(k=0; k<lnz; ++k){
		for(j=0; j<lny; ++j){
			for(i=0; i<lnx; ++i){
				printf("%1.1f ", local_rect[l(i,j,k,lnx,lny)]);
			}
			printf("\n");
		}
		printf("\n\n");
	}

}
/*
if (proc_id == 5){
	
	printf("NW=%d\n",NW);
	MPI_Recv(NW_Recv, lnx*Msp*Msp, MPI_DOUBLE, NW, 99, MPI_COMM_WORLD, &status);
	

	idx[0] = 0; idx[1]=0; idx[2]=isize[2]-Msp;
	dimRbuffer[0]=lnx; dimRbuffer[1]=Msp; dimRbuffer[2]=Msp;
	copyRbuffer(NW_Recv, local_rect, idx, dimRbuffer, isize );

	for(k=0; k<isize[2];++k){
	for(j=0; j<isize[1];++j){
	for(i=0; i<isize[0];++i){
		printf("%1.1f ", local_rect[l(i,j,k,isize[0],isize[1])]);	
	}
	printf("\n");
	}
	printf("\n\n");
	}
	printf("received\n");

}
*/

	Cp3dfft_clean();
	
	free(spread_rect);
	free(local_rect);
	free(N_Recv); free(N_Send); free(S_Recv); free(S_Send); free(W_Recv); free(W_Send);  free(E_Recv); free(E_Send); 
	free(NW_Recv); free(NW_Send); free(SE_Recv); free(SE_Send); free(SW_Recv); free(SW_Send);  free(NE_Recv); free(NE_Send); 


	MPI_Finalize();

	return 0;
}
