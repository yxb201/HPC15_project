all: testfft

testfft: test-fftw3-mpi.c
	mpicc -O3 
