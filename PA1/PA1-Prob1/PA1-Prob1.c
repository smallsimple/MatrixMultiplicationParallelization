#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#define N (2048)
#define threshold (0.000000001)
double A[N][N], B[N][N], RefOut[N][N]; 

int main(){
double rtclock();
void compare();
double clkbegin, clkend;
double t;

int i,j,k,l;
  printf("Matrix Size = %d\n",N);

  for(i=0;i<N;i++)
   for(j=0;j<N;j++)
     {
      A[i][j] = 1.0*(i+0.25*j)/(N+1);
      B[i][j] = 1.0*(i-0.5*j)/(N+1);
     }

  clkbegin = rtclock();
//
// Time the reference baseline version of code 
// Its output is to be compared to the test version to verify correctness
//
  for (i=0; i<N; i++)
   for (j=0; j<N; j++)
    for (k=i; k<N; k++)
      B[i][j] += A[i][k]*B[k][j];
//
// End of reference code
//
  clkend = rtclock();
  t = clkend-clkbegin;
  if (B[N/2][N/2]*B[N/2][N/2] < -1000.0) 
   printf("To foil dead-code elimination by compiler: should never get here\n");
  printf("Base Symm-MatMult: %.1f MFLOPS; Time = %.3f sec; \n",
  1.0*N*N*(N+1)/t/1000000,t);


  for(i=0;i<N;i++)
   for(j=0;j<N;j++)
     {
      RefOut[i][j] = B[i][j];
      A[i][j] = 1.0*(i+0.25*j)/(N+1);
      B[i][j] = 1.0*(i-0.5*j)/(N+1);
     }

  clkbegin = rtclock();


// Test Version      

  for(i=0;i<N;i++)
   for(j=0;j<N;j++)
     {
      A[i][j] = 1.0*(i+0.25*j)/(N+1);
      B[i][j] = 1.0*(i-0.5*j)/(N+1);
     }

  clkbegin = rtclock(); 
//
// Test version of code; initially just contains a copy of base code
// To be modified by you to improve performance
//
  for (i=0; i<N; i++)
   for (j=0; j<N; j++)
    for (k=i; k<N; k++)
      B[i][j] += A[i][k]*B[k][j];
//
  clkend = rtclock();
  t = clkend-clkbegin;
  if (B[N/2][N/2]*B[N/2][N/2] < -1000.0) 
   printf("To foil dead-code elimination by compiler: should never get here\n");
  printf("Test Symm-MatMult: %.1f MFLOPS; Time = %.3f sec; \n",
  1.0*N*N*(N+1)/t/1000000,t);

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
     this_diff = RefOut[i][j]-B[i][j];
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

