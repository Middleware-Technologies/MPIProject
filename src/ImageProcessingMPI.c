/*
 Name        : ImageProcessingMPI.c
 Author      : Manca-Nero-Taina
 Description : Use MPI (+OpenMP) to implement a parallel version of the gamma correction image enhancement 
               algorithm for grayscale images. Use the image format you prefer ("pgm - portable graymap" 
               is an easy to parse grayscale format). The code has to be demonstrated using at least two 
               physical machines connected in a LAN (with or without additional virtual machines to emulate 
               a larger cluster). 

##COMPILAZIONE DA TERMINALE
> mpicc src/ImageProcessingMPI.c -o out/gammaCorrection -lnetpbm -lm -fopenmp


##ESECUZIONE MULTITHREADING DA TERMINALE
> mpirun -np NUM_THREADS out/gammaCorrection IMG_PATH GAMMA_VALUE
E.G.: mpirun -np 2 out/gammaCorrection img/galaxy.ascii.pgm 2


##ESECUZIONE IN CLUSTER
PATH/mpirun -host LIST OUT
/home/middleware/.openmpi/bin/mpirun --host 192.168.0.100,192.168.0.102 out/OUT

*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pgm.h>
#include <math.h>
#include <omp.h>

int main (int argc, char **argv) {

	int i,j;                            /* useful indexes */

	int gamma=2;                        /* correction parameter */
	gray max;                           /* maximum value */

	char *fileName="img/imm.pgm";       /* input file path */
	char *newFileName="out/imm.pgm";    /* output file path */


	gray **image;                       /* to get the output from netpgm's functions */
	int *imageArray;                    /* 2D array containing the image */
	int *procRow;                       /* the data that each processor will process */
	int n_rows, n_cols;                 /* image dimensions */

	int my_rank;                        /* rank of process */
	int n_procs;                        /* number of processes */

	int *sendcounts, *displs;           /* distribution of the numbers among processes */

	// Starting MPI and retrieving some parameters
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	sendcounts = malloc(n_procs*sizeof(int));
	displs = malloc(n_procs*sizeof(int));

	if (my_rank == 0) {

		// Getting the arguments 
		if ( argc == 3 ) {
			fileName = argv[1];
			gamma = atoi(argv[2]);
		}

		// Init netpgm library
		pgm_init(&argc, argv);

		// Get the data from the file
		FILE *f = fopen(fileName,"r");
		image = pgm_readpgm(f, &n_cols, &n_rows, &max);

		// Allocating the data to a contiguous array
		// TODO: maybe this can be avoided but I wasn't able to
		imageArray = malloc(sizeof *imageArray * n_rows * n_cols);
		for (i=0; i<n_rows; i++) {
			for (j=0; j<n_cols; j++){
				imageArray[j*n_rows+i] = image[i][j];
			}
		}

		// Calculating the workload distribution
		// Rows are not important, the total array will be scattered equally...
		int div = (n_rows*n_cols)/n_procs;
		for(i=0; i<n_procs; i++) {
			sendcounts[i] = div;
			displs[i] = i*div;
		}

		// ...except for the process with higher rank, which takes the remains of the division
		sendcounts[n_procs-1] += (n_rows*n_cols) % n_procs;
	}

	// Broadcasting useful values from the root
	MPI_Bcast(&max, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&gamma, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(sendcounts, n_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs, n_procs, MPI_INT, 0, MPI_COMM_WORLD);

	// ... making sure everyone has them before starting
	MPI_Barrier(MPI_COMM_WORLD);

	procRow = malloc(sendcounts[my_rank]*sizeof(int));

	MPI_Scatterv(	imageArray, sendcounts, displs, MPI_INT,
					procRow, sendcounts[my_rank], MPI_INT,
					0, MPI_COMM_WORLD
	);
    
	#pragma omp parallel for
	for(i=0; i<sendcounts[my_rank]; i++) {
		// Applying the correction
		double base = (double)procRow[i]/(double)max;
		procRow[i] = max*pow(base, gamma);
	}

	MPI_Gatherv(    procRow, sendcounts[my_rank], MPI_INT, 
					imageArray, sendcounts, displs, MPI_INT,
					0, MPI_COMM_WORLD
		);

	// Let's wait everyone before writing results
	MPI_Barrier(MPI_COMM_WORLD);

	if (my_rank == 0) {
		// Copy everything back
		for (i=0; i<n_rows; i++) {
			for (j=0; j<n_cols; j++){
				image[i][j] = imageArray[j*n_rows+i];
			}
		}

		// Create the resulting image
		FILE *fout=fopen(newFileName,"w");
		pgm_writepgm(fout, image, n_cols, n_rows, max, 1);

		// Free space
		pgm_freearray(image, n_rows);
		free(imageArray);
	}

	// Release memory
	free(procRow);

	// That's it folks!
	MPI_Finalize();

	return 0;
}