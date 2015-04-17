/*
 Name        : ImageProcessingMPI.c
 Author      : Manca-Nero-Taina
 Description : Use MPI (+OpenMP) to implement a parallel version of the gamma correction image enhancement 
               algorithm for grayscale images. Use the image format you prefer ("pgm - portable graymap" 
               is an easy to parse grayscale format). The code has to be demonstrated using at least two 
               physical machines connected in a LAN (with or without additional virtual machines to emulate 
               a larger cluster). 

##COMPILAZIONE DA TERMINALE
mpicc src/ImageProcessingMPI.c -o out/OUT -I /usr/local/netpbm/include/ -L /usr/local/netpbm/lib -lnetpbm -lm

##ESECUZIONE MULTITHREADING DA TERMINALE
PATH/mpirun -np NUM OUT
 */

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <pgm.h>


int main(int argc, char *argv[])
{
	int my_rank;	    /* rank of process */
	int num_procs;	    /* number of processes */

	/* start up MPI */
	MPI_Init(&argc, &argv);
	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 


	printf("Process: %d - Total Process: %d\n",my_rank,num_procs);


	if (my_rank == 0) //IF I'M MASTER THREAD
	{
	  //GET THE MATRIX
	  //2d array that contains image's data
	  gray **image;

	  //Maximum value of our input image, probably
	  gray max;

	  //Num cols,row,matrix's indices
	  int cols, rows,y,x;

	  //Initialize libpgm
	  pgm_init(&argc, argv);
   
	  //Read the image
	  FILE *f=fopen("imm.pgm","r");
	  image = pgm_readpgm(f, &cols, &rows, &max);
	  int imageMatrix[rows][cols];

	  for (y=0; y<rows; y++)
	    {
	      for (x=0; x<cols; x++)
		{
		  imageMatrix[y][x] = image[y][x];
		  if( imageMatrix[y][x]<10)
		    printf("%d  ",imageMatrix[y][x]);
		  else
		    printf("%d ",imageMatrix[y][x]);
		}
	      printf("\n");
	    }
	  pgm_freearray(image, rows);
	
	  //1 - SEND PORTION OF THE MATRIX TO OTHER THREADS (WITH OPENMP)

	  //2 - PROCESS IMAGE WITH OPENMP FOR INCREMENT PERFORMANCE

	  //3 - ATTENDS PART OF IMAGE FROM OTHER THREAD

	  //4 - CREATE FINAL IMAGE
	  // FILE *fout=fopen("imm2.pgm","w");
	  //pgm_writepgm(fout, image, cols, rows, max, 1);



	}
	else //IF I'M ANOTHER THREAD
	{
		//1 - READ MY PORTION OF IMAGE/RECEIVE PORTION FROM MASTER
		//1 - READ MY PORTION OF IMAGE
		//2 - PROCESS IMAGE WITH OPENMPFOR INCREMENT PERFORMANCE
		//3 - SEND FINAL PART OF IMAGE TO MASTER THREAD
	}

	/* shut down MPI */
	MPI_Finalize(); 
	return 0;
}


