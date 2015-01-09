#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include "../MyMPI.h"
//#define nx 4000
#define nx 10
//#define ny 4000
#define ny 5 
//#define tmax 100
#define tmax 10
#define coeff1 0.5
#define coeff2 0.7
//typedef 
//#define MPI_TYPE MPI_DOUBLE 

double ex[nx][ny];
double ey[nx][ny];
double hz[nx][ny];
double g_hz[nx][ny];
void check();
void dump();
void mpi_check(double **,int ,int );

int main(int argc, char *argv[]){
  double rtclock();
				
  int tt,i,j,nt;
  double clkbegin, clkend;
  double t, maxdiff;

	int id;              // process ID number
	int p;               // number of processes
	MPI_Init (&argc, &argv); 
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);
	MPI_Barrier(MPI_COMM_WORLD);
 
// Initialize arrays
	if(id==p-1){
  	for (i=0; i<nx; i++){
    	for (j=0; j<ny; j++) {
      	ex[i][j] = sin(i)*(1-sin(j));
			  //printf("ex[%d][%d]=%f\t",i,j,ex[i][j]);	
    	}
			//printf("\n");	
  	}
		//printf("\n");
  	for (i=0; i<nx; i++) {
    	for (j=0; j<ny; j++) {
      	ey[i][j] = cos(i)*(1-cos(j));
				//printf("ey[%d][%d]=%f\t",i,j,ey[i][j]);
    	}
			//printf("\n");
  	}
		printf("\n");	
  	for (i=0; i<nx; i++) {
    	for (j=0; j<ny; j++) {
      	hz[i][j] = sin(i)*(1-cos(j));
			//	printf("hz[%d][%d]=%f\t",i,j,hz[i][j]);
    	}
			//printf("\n");	
  	}
		//printf("\n");
		printf("seq:\n");
  	clkbegin = rtclock();
  	for (tt=0; tt<tmax; tt++){
  	  for (j=0; j<ny; j++)
    	  ey[0][j] = tt;
		
     	for (j=0; j<ny; j++) {
      	for (i=1; i<nx; i++){ 
				//	printf("bs1:ey[%d][%d]=%f\t",i,j,ey[i][j]);
       		ey[i][j] = ey[i][j] - coeff1*(hz[i][j]-hz[i-1][j]);
				//	printf("as1:ey[%d][%d]=%f\t",i,j,ey[i][j]);
				}
				//printf("\n");
			}
     	for (j=1; j<ny; j++){ 
      	for (i=0; i<nx; i++){ 
       		ex[i][j] = ex[i][j] - coeff1*(hz[i][j]-hz[i][j-1]);
					//printf("s2:ex[%d][%d]=%f\t",i,j,ex[i][j]);
				}
				//printf("\n");
			}
     	for (j=0; j<ny-1; j++){ 
      	for (i=0; i<nx-1; i++){ 
       		hz[i][j] =  hz[i][j] -
       	  	coeff2*(ex[i][j+1]-ex[i][j]+ey[i+1][j]-ey[i][j]);
					//printf("s3:hz[%d][%d]=%f\t",i,j,hz[i][j]);
				}
				//printf("\n");
			}
  	}
  	clkend = rtclock();
  	t = clkend-clkbegin;
  	printf ("Sequential GFLOPS: %.1f, Time: %.1f\n", 10.0*nx*ny*tmax/t/1e9,t);
  	check();

		for (i=0; i<nx; i++){
			for (j=0; j<ny; j++)
				printf("hz[%d][%d]=%f\t",i,j,hz[i][j]);		
			printf("\n",tt);
		}
	}
	
	if (id==p-1){
		printf("seq end!\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	int nxp;
	int m = nx%p;
	MPI_Status status;
// mpi initialize 
	if (m==0)
		nxp = nx/p;
	else
		{
			if (id<m) nxp = nx/p+1;
			else 	nxp = nx/p;
		}
	printf("id=%d  nxp=%d\n",id,nxp);
  double exp[nxp+2][ny];
	double eyp[nxp+2][ny];
	double hzp[nxp+2][ny];
	int offset;
	if (id>=m) offset=m;
	else offset=0;
  for (i=1; i<nxp+1; i++){
 		for (j=0; j<ny; j++) {
   		exp[i][j] = sin(i-1+id*nxp+offset)*(1-sin(j));
		//	printf("%d exp[%d][%d]=%f\t",id,i-1+id*nxp+offset,j,exp[i][j]);
 		}
		//printf("\n");
 	}
	//printf("\n");
 	for (i=1; i<nxp+1; i++) {
   	for (j=0; j<ny; j++) {
   		eyp[i][j] = cos(i-1+id*nxp+offset)*(1-cos(j));
			//printf("pt%d eyp[%d][%d]=%f\t",id,i-1+id*nxp+offset,j,eyp[i][j]);
   	}
		//printf("\n");
 	}
	//printf("\n");
 	for (i=1; i<nxp+1; i++) {
 		for (j=0; j<ny; j++) {
   		hzp[i][j] = sin(i-1+id*nxp+offset)*(1-cos(j));
		//	printf("%d hzp[%d][%d]=%f\t",id, i-1+id*nxp+offset,j,hzp[i][j]);
   	}
		//printf("\n");
 	}
	//printf("\n");
	MPI_Barrier(MPI_COMM_WORLD);
// calculate
	int is, ie;
	if (id==0) is=2;
	else is=1;
	if (id==p-1) ie=nxp-1;
	else ie=nxp;
		MPI_Barrier(MPI_COMM_WORLD);
		if(id<p-1) MPI_Send(&hzp[nxp],ny,MPI_DOUBLE,id+1,10,MPI_COMM_WORLD);
		if(id>0) MPI_Recv(&hzp[0],ny,MPI_DOUBLE,id-1,10,MPI_COMM_WORLD,&status);
		if(id>0) MPI_Send(&eyp[1],ny,MPI_DOUBLE,id-1,20,MPI_COMM_WORLD);
		if(id<p-1) MPI_Recv(&eyp[nxp+1],ny,MPI_DOUBLE,id+1,20,MPI_COMM_WORLD,&status);
		MPI_Barrier(MPI_COMM_WORLD);
/*
 	for (i=1; i<nxp+1; i++) {
   	for (j=0; j<ny; j++) {
   		eyp[i][j] = cos(i-1+id*nxp+offset)*(1-cos(j));
		//	printf("apt%d eyp[%d][%d]=%f\t",id,i-1+id*nxp+offset,j,eyp[i][j]);
   	}
		//printf("\n");
 	}
*/
  for (tt=0; tt<tmax; tt++){
/*		MPI_Barrier(MPI_COMM_WORLD);
		if(id<p-1) MPI_Send(&hzp[nxp],ny,MPI_DOUBLE,id+1,tt,MPI_COMM_WORLD);
		if(id>0) MPI_Recv(&hzp[0],ny,MPI_DOUBLE,id-1,tt,MPI_COMM_WORLD,&status);
		if(id>0) MPI_Send(&eyp[1],ny,MPI_DOUBLE,id-1,tt+tmax,MPI_COMM_WORLD);
		if(id<p-1) MPI_Recv(&eyp[nxp+1],ny,MPI_DOUBLE,id+1,tt+tmax,MPI_COMM_WORLD,&status);
		MPI_Barrier(MPI_COMM_WORLD);
*/
/*
		if (id==p-1)
			for(j=0; j<ny;j++){
				printf("%d hzp[s][%d]=%f\t",id,j,hzp[0][j]);
			}
		printf("\n");
		if (id==0)
			for(j=0; j<ny; j++)
				printf("%d hzp[n][%d]=%f\t",id,j,hzp[nxp][j]);
		printf("\n");
*/
// stop here
//		exit (0);
		if(id ==0)
 	  	for (j=0; j<ny; j++)
   	  	eyp[1][j] = tt;

   	for (i=is; i<=nxp; i++) 
		{
   		for (j=0; j<ny; j++){ 
				//printf("bp1:%d eyp[%d][%d]=%f\t", id, i-1+offset+id*nxp, j, eyp[i][j]);
      	eyp[i][j] = eyp[i][j] - coeff1*(hzp[i][j]-hzp[i-1][j]);
				//printf("ap1:%d eyp[%d][%d]=%f\t", id, i-1+offset+id*nxp, j, eyp[i][j]);
			}
			//printf("\n");
		}
		//printf("\n");
		MPI_Barrier(MPI_COMM_WORLD);
		if(id>0) MPI_Send(&eyp[1],ny,MPI_DOUBLE,id-1,tt+tmax,MPI_COMM_WORLD);
		if(id<p-1) MPI_Recv(&eyp[nxp+1],ny,MPI_DOUBLE,id+1,tt+tmax,MPI_COMM_WORLD,&status);
		MPI_Barrier(MPI_COMM_WORLD);
		
	/*	if (id==p-1){
			printf("send\t");
			for(j=0; j<ny;j++){
				printf("%d eyp[s][%d]=%f\t",id,j,eyp[1][j]);
			}
		}
		printf("\n");
		if (id==0){
			printf("receive\t");
			for(j=0; j<ny; j++){
				printf("%d eyp[r][%d]=%f\t",id,j,eyp[nxp+1][j]);
			}
		}
		printf("\n");
		*/
    for (i=1; i<=nxp; i++){ 
   		for (j=1; j<ny; j++){ 
       	exp[i][j] = exp[i][j] - coeff1*(hzp[i][j]-hzp[i][j-1]);
				//printf("p2:%d exp[%d][%d]=%f\t",id,i-1+offset+id*nxp,j,eyp[i][j]);
			}
			//printf("\n");
		}
    for (i=1; i<=ie; i++){
   		for (j=0; j<ny-1; j++){ 
       	hzp[i][j] =  hzp[i][j] -
       	  coeff2*(exp[i][j+1]-exp[i][j]+eyp[i+1][j]-eyp[i][j]);
				//printf("p3:%d hzp[%d][%d]=%f\t",id,i-1+offset+id*nxp,j,hz[i][j]);
			}
			//printf("\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(id<p-1) MPI_Send(&hzp[nxp],ny,MPI_DOUBLE,id+1,tt,MPI_COMM_WORLD);
		if(id>0) MPI_Recv(&hzp[0],ny,MPI_DOUBLE,id-1,tt,MPI_COMM_WORLD,&status);
		MPI_Barrier(MPI_COMM_WORLD);
 	}
	if(id<p-1) MPI_Send(&hzp[nxp],ny,MPI_DOUBLE,id+1,tt,MPI_COMM_WORLD);
	if(id>0) MPI_Recv(&hzp[0],ny,MPI_DOUBLE,id-1,tt,MPI_COMM_WORLD,&status);

	MPI_Barrier(MPI_COMM_WORLD);

	for (i=1; i<=nxp; i++){
		for (j=0; j<ny; j++)
			printf("%d hzp[%d][%d]=%f\t",id,i-1+id*nxp+offset,j,hzp[i][j]);
		printf("\n");
	}

//	MPI_Barrier(MPI_COMM_WORLD);
// check
	//printf("bobo1, id=%d\n",id);
	MPI_Barrier(MPI_COMM_WORLD);
/*
	for (i=1; i<=nxp; i++){
		for (j=0; j<ny; j++)
			g_hz[i-1+id*nxp+offset][j]=hzp[i][j];
		MPI_Bcast(&g_hz[i-1+id*nxp+offset],ny,MPI_DOUBLE,id,MPI_COMM_WORLD);
	}
  MPI_Barrier(MPI_COMM_WORLD);
	if (id ==0)
		for (i=0; i<nx; i++){
			for (j=0; j<ny; j++)
				printf("g_hz[%d][%d]=%f\t",i,j,g_hz[i][j]); 
			printf("\n");
		}
*/

	int count1,count2;
	count1=0;
	count2=0;
	if (id>0)
		for (i=1; i<=nxp; i++){
			//printf("bobo1.01\n");
			MPI_Send(&hzp[i],ny,MPI_DOUBLE,0,i-1+id*nxp+offset,MPI_COMM_WORLD);
			printf("bobo:send  id=%d, i=%d, n_i=%d\n",id,i,i-1+id*nxp+offset);
			count1++;
		}

	//MPI_Barrier(MPI_COMM_WORLD);	
	//printf("bobo1.1 id=%d\n",id);
	if (p>1 && id==0){
		printf("bobo2,nxp=%d, nx=%d\n",nxp,nx);
		int idr;
		for (i=nxp; i<nx; i++){
			if (i<m) idr=i/(nxp+1);
			else idr= (i-m)/nxp + m;
			printf("bobo3.9:b recv i=%d idr=%d\n",i,idr);
			MPI_Recv(&g_hz[i],ny,MPI_DOUBLE,idr,i,MPI_COMM_WORLD,&status);
			printf("bobo4:a Recv i=%d, idr=%d\n",i, idr);
			count2++;
		}
		printf("count1=%d  count2=%d\n",count1,count2);
		printf("bobo4.5 \n");
		MPI_Barrier(MPI_COMM_WORLD);
		for (i=0; i< nxp; i++){
			for (j=0; j<ny; j++)
				g_hz[i][j]=hzp[i+1][j];
		}		

		printf("bobo5\n");
		for (i=0; i<nx; i++){
			for (j=0; j<ny; j++)
				printf("g_hz[%d][%d]=%f\t",i,j,g_hz[i][j]);
			printf("\n");
		}
	}


}
	

double rtclock()
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

void check()
{
double maxval, minval, mean, rms;
int i,j;
  minval = maxval = hz[0][0];
  mean = rms = 0.0;
  for (i=0;i<nx;i++)
   for (j=0;j<ny;j++)
    {
     if (hz[i][j] < minval) minval = hz[i][j];
     if (hz[i][j] > maxval) maxval = hz[i][j];
     mean += hz[i][j]/(1.0*nx*ny);
     rms += hz[i][j]*hz[i][j]/(1.0*nx*ny);
    }
  rms = sqrt(rms);
  printf("Minhz= %18.9f; Maxhz = %18.9f; Mean= %18.9f, RMS = %18.9f\n", minval,maxval,mean,rms);
}


void mpi_check(double **hz_g,int nxp_g, int ny_g)
{
	int i,j;
	double maxval, minval, mean, rms;
	double g_maxval, g_minval, g_mean, g_rms;
	
	minval = maxval = hz_g[0][0];
	mean = rms = 0.0;
	for (i=0; i< nxp_g; i++)
		for(j=0; j<ny_g; j++)
		{
			if (hz_g[i][j] < minval) minval = hz_g[i][j];
			if (hz_g[i][j] > maxval) maxval = hz_g[i][j];
			mean += hz_g[i][j]/(1.0*nx*ny);
			rms += hz_g[i][j];
		}



}
