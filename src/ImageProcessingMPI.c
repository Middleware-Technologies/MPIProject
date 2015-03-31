/*
 ============================================================================
 Name        : ImageProcessingMPI.c
 Author      : Manca-Nero-Taina
 Version     :
 Copyright   : 
 Description : Calculate Pi in MPI
 ============================================================================
 */
#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
	int			my_rank;	    	/* rank of process */
	int			num_procs;		    /* number of processes */


	/* start up MPI */
	MPI_Init(&argc, &argv);
	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 


	printf("Process: %d - Total Process: %d\n",my_rank,num_procs);


	if (my_rank == 0) //IF I'M MASTER THREAD
	{


		//(SEND PART OF IMAGE TO OTHER THREAD)
		//1 - READ MY PORTION OF IMAGE
		//2 - PROCESS IMAGE WITH OPENMPFOR INCREMENT PERFORMANCE
		//3 - ATTENDS PART OF IMAGE FROM OTHER THREAD
		//4 - CREATE FINAL IMAGE


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
