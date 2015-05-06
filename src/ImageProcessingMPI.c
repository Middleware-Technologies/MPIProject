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
mpirun -np NUM OUT

##ESECUZIONE IN CLUSTER
PATH/mpirun -host LIST OUT
/home/middleware/.openmpi/bin/mpirun --host 192.168.0.100,192.168.0.102 out/OUT

Read: http://stackoverflow.com/questions/9269399/sending-blocks-of-2d-array-in-c-using-mpi
*/

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <pgm.h>
#include <math.h>


int applyCorrection(int input, double gamma, int depth)
{
	double base=(double)input/(double)depth;
	return (int)(depth* pow(base,1/gamma));
}

int main(int argc, char *argv[]) 
{
	int my_rank; /* rank of process */
	int num_procs; /* number of processes */
       
	/* start up MPI */
	MPI_Init(&argc, &argv);
	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	printf("Process: %d - Total Process: %d\n", my_rank, num_procs);

	if (my_rank == 0) {		// IF I'M MASTER THREAD
	        //Data Structure and utilities
	        gray **imageArray;
		gray max;      
	        unsigned int row, column,x,y;

		//0 - READ FILE'S MATRIX
		pgm_init(&argc, argv);

		// Read the input image and pass
		FILE *f = fopen("img/imm.pgm", "r");
		imageArray = pgm_readpgm(f, &column, &row, &max);
		
		//TEST: Print Matrix
		printf("Thread 0 Read:\n");
		for (x = 0; x < row; x++) 
		{
			for (y = 0; y < column; y++) 
			{				  				 
				if(imageArray[x][y]<10)
					printf("%d  ",imageArray[x][y]);
				else
					printf("%d ",imageArray[x][y]);		  
			}
			printf("\n");
		}
		printf("\n");

		//1 - SEND PORTION OF THE MATRIX TO OTHER THREADS - TEST WITH 2 PROCESS		
		//How many rows we need to send to every process?
		int numRow=3;
		
		int index;  
		for(index=1;index<num_procs;index++) 
		{
			//First: send packet's dimension {numRows,numCols,depth}
			int param[3]={numRow,column,(int)max};
			MPI_Send(param,numRow,MPI_INT,index,1,MPI_COMM_WORLD);

			//Send Packet
			MPI_Send(&imageArray[0][0],numRow*column,MPI_INT,index,2,MPI_COMM_WORLD);
		}

		//2 - PROCESS IMAGE WITH OPENMP FOR INCREMENT PERFORMANCE

		//3 - ATTENDS PART OF IMAGE FROM OTHER THREAD AND UPDATE imageArray MATRIX
		
		//4 - CREATE FINAL IMAGE

		// Free space
		pgm_freearray(imageArray, row);
	} 
	else 
	{
		//1 - RECEIVE PORTION FROM MASTER
		MPI_Status stat;
		
		int param[3];
		MPI_Recv(param,3,MPI_INT,0,1,MPI_COMM_WORLD,&stat);

		int packet[param[0]][param[1]];
		MPI_Recv(packet,param[0]*param[1],MPI_INT,0,2,MPI_COMM_WORLD,&stat);

		//2 - PROCESS IMAGE WITH OPENMPFOR INCREMENT PERFORMANCE
		int i,j;		
		for(i=0;i<param[0];i++)
		{
			for(j=0;j<param[1];j++)
			{
				if(packet[i][j]<10)
			  		printf("%d  ",packet[i][j]);
				else
			  		printf("%d ",packet[i][j]);
				int val=applyCorrection(packet[i][j],2,param[2]);
				packet[i][j]=val;
			}
			printf("\n");
		}
		

		//3 - SEND FINAL PART OF IMAGE TO MASTER THREAD
		//...
	}

	/* shut down MPI */
	MPI_Finalize();

	return 0;
}
