#include <omp.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#define N (16)
#define threshold (0.000000001)
void compare();
double A[N][N][N],C[N][N],RefOut[N][N];

int main(){
  int i,j,k,l;
  int nthr;

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	A[i][j][k]=2i-j-k;

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      C[i][j]=i-j;
  //  C=1;

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	C[i][k]+=A[i][j][k];

    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            RefOut[i][j] = C[i][j];


    //  printf("orignal C\n");
    //for(i=0;i<N;i++)
    //printf("%4f ",C[i]);

  
  nthr=4;
  
  omp_set_num_threads(nthr);

  //  for(i=0;i<N;i++)
  //for(j=0;j<N;j++)
  //  A[i]=i;
  //  for(i=0;i<N;i++)
  //C[i]=i;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      C[i][j]=i-j;
  
#pragma omp parallel private(i,j,k)
  //#pragma omp for
  //#pragma omp critical
#pragma omp for
    for(i=0;i<N;i++)
      //printf("i=%d, threadID=%d\n",i,omp_get_thread_num());
      //#pragma omp for
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	C[i][k]+=A[i][j][k];
    
#pragma omp critical 
  //  {
  // printf("\nparallel C\n");
  // for(i=0;i<N;i++)
  //   printf("%4f ",C[i]);
 
  //  }
  
  compare();
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

