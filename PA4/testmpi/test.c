#include <stdio.h>
#include <math.h>
#include "mpi.h"
int main (int argc, char *argv[]){
	int n;
	int *n_all;
	int p, id;
	double g_data[4];
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	n=0;
	if (id==0)
		n=10;
	printf("1: id=%d, n=%d\n", id, n);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (id==0) printf("after MPI_Bcast\n");
	printf("2: id=%d, n=%d\n", id, n);
	MPI_Barrier(MPI_COMM_WORLD);
	if (id==0) printf("after MPI_Barrier\n");
	MPI_Barrier(MPI_COMM_WORLD);
	n=n+id;
	printf("3: id=%d, n=%d\n", id,n);
	int size,i,j;
	size=4;
	double data[size][size];
	for (i=0; i<size; i++)
		for (j=0; j<size; j++)
		 	data[i][j]=i*10+j+id*100;
	//printf("n = n + id\n");
	if (id ==1) printf("1. id=%d\n",id);
		for (i=0; i<size; i++){
			for (j=0; j<size; j++)
				printf("s%d data[%d][%d]=%f\t",id,i,j,data[i][j]);
			printf("\n");
		}
	MPI_Barrier(MPI_COMM_WORLD);
	//if (id==0) MPI_Send(&data[2],4,dMPI_DOUBLE,1,0,MPI_COMM_WORLD);
	//if (id==1) MPI_Recv(&data[2],4,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&data[1],&g_data,4,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	//MPI_Bcast (&data[1],4,MPI_DOUBLE,1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (id==0)
	{	
		for (i = 0; i<size; i++)
			printf("g_data[%d]=%f\t",i, g_data[i]);
	}
	printf("\n");
	MPI_Finalize();	
}
