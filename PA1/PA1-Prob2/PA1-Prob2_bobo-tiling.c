#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#define N (128)
#define threshold (0.000000001)
double A[N][N][N], B[N][N][N], C[N][N],  RefOut[N][N]; 

int main(){
double rtclock();
void compare();
double clkbegin, clkend;
double t;

int i,j,k,l;
 int it,jt,kt,lt;
  printf("Tensor Size = %d\n",N);

  for(i=0;i<N;i++)
   for(j=0;j<N;j++)
    for (k=0;k<N;k++)
     {
      A[i][j][k] = 1.0*(i+0.25*j+0.5*k)/(N+1);
      B[i][j][k] = 1.0*(i-0.5*j+0.25*k)/(N+1);
     }

  for(i=0;i<N;i++)
   for(j=0;j<N;j++)
    C[i][j] = 0;

  clkbegin = rtclock();
//
// Time the reference baseline version of code 
// Its output is to be compared to the test version to verify correctness
//
  for (i=0; i<N; i++)
   for (j=0; j<N; j++)
    for (k=0; k<N; k++)
     for (l=0; l<N; l++)
       C[i][j] += A[l][i][k]*B[k][j][l];
//
// End of reference code
//
  clkend = rtclock();
  t = clkend-clkbegin;
  if (C[N/2][N/2]*C[N/2][N/2] < -1000.0) 
   printf("To foil dead-code elimination by compiler: should never get here\n");
  printf("Base-TensorMult: %.1f MFLOPS; Time = %.3f sec; \n",
  2.0*N*N*N*N/t/1000000,t);

  for(i=0;i<N;i++)
   for(j=0;j<N;j++)
    RefOut[i][j] = C[i][j];



//
// Initialization for test version of code
//
  for(i=0;i<N;i++)
   for(j=0;j<N;j++)
     C[i][j] = 0;

  clkbegin = rtclock(); 
//
// Test version of code; initially just contains a copy of base code
// To be modified by you to improve performance
//
  for (lt=0;lt<N;lt=lt+N/8)
    for (kt=0;kt<N;kt=kt+N/8)
      for (i=0; i<N; i++)
	for (j=0; j<N; j++)
	  for (k=kt; k<kt+N/8; k++)
	    for (l=lt; l<lt+N/8; l++)
	    {
	      C[i][j] += A[l][i][k]*B[k][j][l];
	    }
//
  clkend = rtclock();
  t = clkend-clkbegin;
  if (C[N/2][N/2]*C[N/2][N/2] < -1000.0) 
   printf("To foil dead-code elimination by compiler: should never get here\n");
  printf("Test-TensorMult: %.1f MFLOPS; Time = %.3f sec; \n",
  2.0*N*N*N*N/t/1000000,t);
//
// Verify correctness by comparing result with reference version's
//
  compare();
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

void compare()
{
double maxdiff,this_diff;
int numdiffs;
int i,j;
  numdiffs = 0;
  maxdiff = 0;
  for (i=0;i<N;i++)
   for (j=0;j<N;j++)
    {
     this_diff = RefOut[i][j]-C[i][j];
     if (this_diff < 0) this_diff = -1.0*this_diff;
     if (this_diff>threshold)
      { numdiffs++;
        if (this_diff > maxdiff) maxdiff=this_diff;
      }
    }
   if (numdiffs > 0)
      printf("%d Diffs found over threshold %f; Max Diff = %f\n",
               numdiffs,threshold,maxdiff);
   else
      printf("No differences found between base and test versions\n");
}

