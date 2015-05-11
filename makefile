all: randc randf

# cims
randc: driver_rand.c
	mpicc -I/usr/local/pkg/p3dfft/2.7.1/include -c driver_rand.c
	mpifort -o randc driver_rand.o -L/usr/local/pkg/p3dfft/2.7.1/lib -lp3dfft -lfftw3
randf: driver_rand.F90
	mpifort -I/usr/local/pkg/p3dfft/2.7.1/include -c driver_rand.F90
	mpifort -o randf driver_rand.o -L/usr/local/pkg/p3dfft/2.7.1/lib -lp3dfft -lfftw3
sinef: driver_sine.F90
	mpifort -I/usr/local/pkg/p3dfft/2.7.1/include -c driver_sine.F90
	mpifort -o sinef driver_rand.o -L/usr/local/pkg/p3dfft/2.7.1/lib -lp3dfft -lfftw3

clean:
	rm -f randc randf sinef
