/*
 Name        : ImageProcessingMPI.c
 Author      : Manca-Nero-Taina
 Description : Use MPI (+OpenMP) to implement a parallel version of the gamma correction image enhancement 
               algorithm for grayscale images. Use the image format you prefer ("pgm - portable graymap" 
               is an easy to parse grayscale format). The code has to be demonstrated using at least two 
               physical machines connected in a LAN (with or without additional virtual machines to emulate 
               a larger cluster). 

##COMPILAZIONE DA TERMINALE
mpicc src/ImageProcessingMPI.c -o out/OUT -I /usr/local/netpbm/include/ -L /usr/local/netpbm/lib -lnetpbm -lm -fopenmp

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
#include <omp.h>

int applyCorrection(int input, double gamma, int depth)
{
	double base=(double)input/(double)depth;
	return (int)(depth* pow(base,1/gamma));
}

int main(int argc, char *argv[]) 
{
	int gammaParameter=2;
	char *fileName="img/imm2.pgm";
	char *newFileName="img/imm2New.pgm";


	int my_rank;        /* rank of process */
	int num_procs;      /* number of processes */
       
	
	MPI_Init(&argc, &argv);          	        /* start up MPI */	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);	/* find out process rank */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);      /* find out number of processes */

	
	if (my_rank == 0) {		// IF I'M MASTER THREAD
	        //Data Structure and utilities
	        gray **imageArray;
		gray max;      
	        int row, column,x,y;


		//0 - READ FILE'S MATRIX
		pgm_init(&argc, argv);

		// Read the input image and pass
		FILE *f = fopen(fileName,"r");
		imageArray = pgm_readpgm(f, &column, &row, &max);
		
		//TEST: Initial Matrix
		/*
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
		printf("\n \n");
		*/
		
		
		//1 - SEND PORTION OF THE MATRIX TO OTHER THREADS	
		int numRow=round((double)row/num_procs);
		
		int index;  
		for(index=1;index<num_procs;index++) 
		{
			//First: send packet's dimension {numRows,numCols,depth}
			int param[3]={numRow,column,(int)max};
			MPI_Send(param,3,MPI_INT,index,1,MPI_COMM_WORLD);

			//Send Packet
			int startRow=(index-1)*(numRow);
			MPI_Send(&imageArray[startRow][0],numRow*column,MPI_INT,index,2,MPI_COMM_WORLD);
		}


		//2,3 - ATTENDS PART OF IMAGE FROM OTHER THREAD AND UPDATE imageArray MATRIX		
		#pragma omp parallel num_threads(num_procs)   /* Add -fopenmp in mpicc command */
		{
			int my_num=omp_get_thread_num();
			if(my_num==0)    /* Correct remaining rows */
			{
				int startRow=(index-1)*(numRow);
				for (x = startRow; x < row; x++) 
				{
					for (y = 0; y < column; y++) 
					{				  				 
						int val=applyCorrection(imageArray[x][y],gammaParameter,(int)max);
						imageArray[x][y]=val;	  
					}
				}
			}
			else		/* Receive correct portions */
			{
				MPI_Status stat;
				int startRow=(my_num-1)*(numRow);
				MPI_Recv(&imageArray[startRow][0],numRow*column,MPI_INT,my_num,my_num,MPI_COMM_WORLD,&stat);
			}
		}
		
		//TEST: Final Matrix
		/*
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
		*/
		
		//4 - CREATE FINAL IMAGE
		FILE *fout=fopen(newFileName,"w");
    		pgm_writepgm(fout, imageArray, column, row, max, 1);
			
		// Free space
		pgm_freearray(imageArray, row);
	} 
	else 
	{
		//1 - RECEIVE PORTION FROM MASTER
		MPI_Status stat;
		
		int param[3];
		MPI_Recv(param,3,MPI_INT,0,1,MPI_COMM_WORLD,&stat);

		int imageArray[param[0]][param[1]];
		MPI_Recv(imageArray,param[0]*param[1],MPI_INT,0,2,MPI_COMM_WORLD,&stat);

		//2 - PROCESS IMAGE WITH OPENMPFOR INCREMENT PERFORMANCE
		int x,y;	
		for(x=0;x<param[0];x++)
		{
			for(y=0;y<param[1];y++)
			{		
				int val=applyCorrection(imageArray[x][y],gammaParameter,param[2]);
				imageArray[x][y]=val;	
			}
		}
		
		//3 - SEND FINAL PART OF IMAGE TO MASTER THREAD
		MPI_Send(&imageArray,param[0]*param[1],MPI_INT,0,my_rank,MPI_COMM_WORLD);
	}

	/* shut down MPI */
	MPI_Finalize();

	return 0;
}
