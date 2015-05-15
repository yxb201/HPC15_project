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
	int lnx, lny, lnz, lkx, lky, lkz;
	int istart[3], isize[3], iend[3];
	int fstart[3], fsize[3], fend[3];
	double *spread_rect, *local_rect, *output_rect;
	double *E2xl, *E2yl, *E2zl;	
	int Msp; 
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
	unsigned char op_f[3]="fft";

	double V0,V1,V2,V3; 

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Status status;

	// 2x2 processor for now
	dims[0]=2; dims[1]=2;
	P1 = dims[0]; P2 = dims[1];	

	// number of sources in each processor
	n_src = 8*8*16;

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
	NW    = mod(k+1,P2)*P1 + mod(j-1,P1);

	// assume for now 
	Mr  = 32; 
	Msp = 12;
	R   = 2;
	M  = Mr/R; 
	tau = (1.*Msp) / (M*M); 

	h = L / Mr;

	// precompute E3 and E4
	double *E3, *E4; 
	E3 = (double *) malloc( sizeof(double)* (Msp+1) );
	E4 = (double *) malloc( sizeof(double)* (M/2+1) );
	for (i=0; i<=Msp ; ++i){
		E3[i]=exp(-(M_PI*i/Mr)*(M_PI*i/Mr)/tau);
	}
	
	for (i=0; i<=M/2 ; ++i){
		E4[i] = exp(tau*i*i);
	}
/*
if (proc_id == 0){

//	printf("%f\n", E3[0]);

	printf("E3: ");
	for (i=0; i<=Msp ; ++i){
		printf("%f ",E3[i]);
	}
	printf("\n");

	printf("E4: ");
	for (i=0; i<=M/2 ; ++i){
		printf("%f ",E4[i]);
	}
	printf("\n");

}
*/
	E2xl = (double *) malloc( sizeof(double) * 2*Msp );
	E2yl = (double *) malloc( sizeof(double) * 2*Msp );
	E2zl = (double *) malloc( sizeof(double) * 2*Msp );

	// initialize P3DFFT
	Cp3dfft_setup(dims,Mr,Mr,Mr,MPI_Comm_c2f(MPI_COMM_WORLD), Mr,Mr,Mr, 0, memsize);
	
	// set input dimensions	
	conf = 1;
	Cp3dfft_get_dims(istart, iend, isize, conf);

	
	// set output dimensions
	conf = 2;
	Cp3dfft_get_dims(fstart, fend, fsize, conf);
/*
if (proc_id == 0){
	

	printf("istart: %d %d %d \n", istart[0], istart[1], istart[2]);
	printf("iend: %d %d %d \n", iend[0], iend[1], iend[2]);
	printf("isize: %d %d %d \n", isize[0], isize[1], isize[2]);

	printf("\n");

	printf("fstart: %d %d %d \n", fstart[0], fstart[1], fstart[2]);
	printf("fend: %d %d %d \n", fend[0], fend[1], fend[2]);
	printf("fsize: %d %d %d \n", fsize[0], fsize[1], fsize[2]);

}
*/	
	// dimension of local rectangle
	lnx = isize[0];
	lny = isize[1];
	lnz = isize[2];

	lkx = fsize[0];
	lky = fsize[1];
	lkz = fsize[2];
	
	// one source for now
	xj[0]= (istart[0]+iend[0])*h/2; 
	yj[0]= (istart[1]+iend[1])*h/2; 
	zj[0]= (istart[2]+iend[2])*h/2;

/*	
if (proc_id == 0){
	printf("%f, %f, %f\n",xj[0], yj[0], zj[0]);
}
*/

	// rectangle for spreading, dimension: nx x (ny_local + 2Msp) x (nz_local+2Msp)
	spread_rect = (double *) malloc( sizeof(double) * lnx*(lny+2*Msp)*(lnz+2*Msp) );
	for (i=0; i<lnx*(lny+2*Msp)*(lnz+2*Msp); ++i) spread_rect[i] = 0.;

	dimSpreadRect[0]=lnx; dimSpreadRect[1]=lny+2*Msp; dimSpreadRect[2] = lnz+2*Msp;

	// rectangle for local data, dimension: nx x ny_local x nz_local
	local_rect  = (double *) malloc( sizeof(double) * lnx*lny*lnz );	
	for (i=0; i<lnx*lny*lnz; ++i) local_rect[i] = 0.;

	// set dimension of output data
	output_rect = (double *) malloc( sizeof(double) * fsize[0]*fsize[1]*fsize[2]*2 );


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

	// generate sources
	for (k=0; k<8; k++){
		for(j=0; j<8; j++){
			for (i=0; i<16; i++){
				xj[l(i,j,k,16,8)] = i*(2*M_PI/16);
				yj[l(i,j,k,16,8)] = (j+istart[1]-1)*(2*M_PI/16);
				zj[l(i,j,k,16,8)] = (k+istart[2]-1)*(2*M_PI/16);	
			}			
		}
	}



	// need to print sources to check


	// for each source
	double mxh, myh, mzh;
	double piMtau = M_PI / (Mr * tau);
	for(s=0; s < n_src ; ++s){
		
		// find the closest grid point (in the whole domain)
		mx = (int) ( xj[s]/h );
		my = (int) ( yj[s]/h );
		mz = (int) ( zj[s]/h );
/*		
if (proc_id == 0){
		printf("center: %d %d %d \n", mx, my, mz);
}*/

		mxh = mx*h; myh = my*h; mzh = mz*h;


/*
if (proc_id == 0){
		printf("center: %.16f %.16f %.16f \n", mxh, myh, mzh);
}
		
		// closest grid point (in spreading rect with halo cells )
		smx= mx - (istart[0]-1);
		smy= my - (istart[1]-1) + Msp;
		smz= mz - (istart[2]-1) + Msp;
		
		diffx = xj[s] - mxh;
		diffy = yj[s] - myh;
		diffz = zj[s] - mzh;
		E1 = exp( -(diffx*diffx+diffy*diffy+diffz*diffz)/(4*tau) );

if (proc_id == 0){
		printf("E1 = %.16f \n ", E1);
}

		E2x = exp( piMtau * diffx  );
		E2y = exp( piMtau * diffy  );
		E2z = exp( piMtau * diffz  );


if (proc_id == 0){
		printf("E2x = %.16f, E2y = %.16f, E2z = %.16f \n", E2x, E2y, E2z);
}



		E2xl[Msp-1]=1.; E2yl[Msp-1]=1.; E2zl[Msp-1]=1.;
		for (l1 = 1; l1<=Msp; ++l1 ){
			E2xl[l1+(Msp-1)] = E2xl[(l1-1)+(Msp-1)] * E2x;
			E2yl[l1+(Msp-1)] = E2yl[(l1-1)+(Msp-1)] * E2y; 
 			E2zl[l1+(Msp-1)] = E2zl[(l1-1)+(Msp-1)] * E2z; 
		}

		for (l1 = -1; l1>=-Msp+1; --l1){
			E2xl[l1+(Msp-1)] = E2xl[(l1+1)+(Msp-1)] / E2x; 
			E2yl[l1+(Msp-1)] = E2yl[(l1+1)+(Msp-1)] / E2y; 
			E2zl[l1+(Msp-1)] = E2zl[(l1+1)+(Msp-1)] / E2z; 
		
		}

if (proc_id == 0){
		printf("E2xl: ");
		for (l1 = 0; l1 < 2*Msp; ++l1){
			printf("%.16f ", E2xl[l1]);
		}
		printf("\n");


		printf("E2yl: ");
		for (l1 = 0; l1 < 2*Msp; ++l1){
			printf("%.16f ", E2yl[l1]);
		}
		printf("\n");


		printf("E2zl: ");
		for (l1 = 0; l1 < 2*Msp; ++l1){
			printf("%.16f ", E2zl[l1]);
		}
		printf("\n");
}
*/

		// build the spreading rectangle
		V0 = 1. * E1;
		for (k=-Msp+1; k<=Msp; ++k){

			V1 = V0 * E2zl[k+(Msp-1)] * E3[abs(k)];
			for (j=-Msp+1; j<=Msp; ++j){
				
				V2 = V1 * E2yl[j+(Msp-1)] * E3[abs(j)];
				for (i=-Msp+1; i<=Msp; ++i){ 
					
					V3 = V2 * E2xl[i+(Msp-1)] * E3[abs(i)];
					spread_rect[l(mod(i+smx,Mr),j+smy,k+smz,lnx,lny+2*Msp)] += V3;
				}
			}
		}


	}// end of looping sources 


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

/*
if (proc_id == 0){ 

	for(k=0; k<lnz+2*Msp; ++k){
		for(j=0; j<lny+2*Msp; ++j){
			for(i=0; i<lnx; ++i){
				printf("%1.12f ", spread_rect[l(i,j,k,lnx,lny+2*Msp)]);
			}
			printf("\n");
		}
		printf("\n\n");
	}



	for(k=0; k<lnz; ++k){
		for(j=0; j<lny; ++j){
			for(i=0; i<lnx; ++i){
				printf("%1.12f ", local_rect[l(i,j,k,lnx,lny)]);
			}
			printf("\n");
		}
		printf("\n\n");
	}

}*/

	// step 2: take FFT on local_rect
	MPI_Barrier(MPI_COMM_WORLD);
	Cp3dfft_ftran_r2c(local_rect, output_rect, op_f);

if (proc_id == 1){
	
	for(k = 0; k<lkz; ++k){
		for(j=0; j<lky; ++j){
			for(i=0; i<2*lkx; ++i){
				printf("%1.12f ", output_rect[l(i,j,k,lkx*2,lky)]);
			}
			printf("\n");
		}
		printf("\n \n");
	}

}
	// step 3: Deconvolution

	Cp3dfft_clean();
	
	free(xj); free(yj); free(zj);
	free(E2xl); free(E2yl); free(E2zl);
	free(E3); free(E4);	
	free(spread_rect);
	free(local_rect);
	free(output_rect);
	free(N_Recv); free(N_Send); free(S_Recv); free(S_Send); free(W_Recv); free(W_Send);  free(E_Recv); free(E_Send); 
	free(NW_Recv); free(NW_Send); free(SE_Recv); free(SE_Send); free(SW_Recv); free(SW_Send);  free(NE_Recv); free(NE_Send); 


	MPI_Finalize();

	return 0;
}
