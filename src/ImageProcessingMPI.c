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

 Read: http://stackoverflow.com/questions/9269399/sending-blocks-of-2d-array-in-c-using-mpi
 */

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <pam.h>

int main(int argc, char *argv[]) 
{
	int my_rank; /* rank of process */
	int num_procs; /* number of processes */

	// The two pam, one for input and one for output,structures
	// containing infos about the image
	struct pam inpam, outpam;

	// The 3D array containing the values
	tuple ** imageArray;

	// Indexes to browse the array
	unsigned int row, column, plane;

	/* start up MPI */
	MPI_Init(&argc, &argv);
	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	printf("Process: %d - Total Process: %d\n", my_rank, num_procs);

	if (my_rank == 0) {		// IF I'M MASTER THREAD
		//0 - READ FILE'S MATRIX
		pgm_init(&argc, argv);

		// Read the input image and pass
		inpam.file = fopen("img/imm.pgm", "r");
		imageArray = pnm_readpam(inpam.file, &inpam, PAM_STRUCT_SIZE(tuple_type));
		
		int imageIntArray[inpam.height][inpam.width];
		
		//TEST: Print Matrix
		for (row = 0; row < inpam.height; row++) {
			for (column = 0; column < inpam.width; column++) {
				for (plane = 0; plane < inpam.depth; ++plane) {
				  int val=(int)imageArray[row][column][plane];
				  imageIntArray[row][column]=val;
				  printf("%d ",val);
				}
			}
			printf("\n");
		}
		printf("\n");printf("\n");
	       

		//1 - SEND PORTION OF THE MATRIX TO OTHER THREADS  - TEST WITH 2 PROCESS
		int index;
		int numRows;   

		if(inpam.height%2==0)
		  numRows=inpam.height/2;
		else
		  numRows=(inpam.height-1)/2;
		  
		for(index=1;index<num_procs;index++)
		{
		  //Create a packet to send to index thread
		  int rows=(inpam.height-numRows);
		  int columns=(inpam.width);
		  int param[2]={rows,columns};
		  
		  //First: send packet's dimension
		  MPI_Send(param,2,MPI_INT,index,1,MPI_COMM_WORLD);
		  //Send Packet
		  MPI_Send(&imageIntArray[numRows][0],dimension,MPI_INT,index,2,MPI_COMM_WORLD);
 
		}
		

		//2 - PROCESS IMAGE WITH OPENMP FOR INCREMENT PERFORMANCE

		//3 - ATTENDS PART OF IMAGE FROM OTHER THREAD
		
		//4 - CREATE FINAL IMAGE
		outpam = inpam;
		outpam.file = fopen("img/imm2.pgm", "w");
		outpam.plainformat = 1; // Force plain, so that we, useless meatbags, can read what is being written

		// Write all the resulting values into the file
		pnm_writepam(&outpam, imageArray);

		// Free space
		pnm_freepamarray(imageArray, &inpam);

	} 
	else 
	{
	  //1 - RECEIVE PORTION FROM MASTER
	  MPI_Status stat;
	  int param[2];
	  MPI_Recv(param,2,MPI_INT,0,1,MPI_COMM_WORLD,&stat);

	  int packet[param[0]][param[1]];
	  MPI_Recv(packet,param[0]*param[1],MPI_INT,0,2,MPI_COMM_WORLD,&stat);

	  int i,j;
	  printf("RECEIVED\n");
	  for(i=0;i<param[0];i++)
	  {
	    for(j=0;j<param[1];j++)
	      {
	      printf("%d ",packet[i][j]);
	      }
	    printf("\n");
	  }
	  printf("\n");printf("\n");printf("\n");printf("\n");


	  //2 - PROCESS IMAGE WITH OPENMPFOR INCREMENT PERFORMANCE
	  //3 - SEND FINAL PART OF IMAGE TO MASTER THREAD
	}


	/* shut down MPI */
	MPI_Finalize();

	return 0;
}
