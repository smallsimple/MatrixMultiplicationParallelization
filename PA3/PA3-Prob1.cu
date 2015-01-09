#include <cuda.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>

#define N (2048)
// Use N=2048 for your testing
//#define N (256) 
//#define N (256)
#define threshold (0.000000001)
double rtclock();
void compare();
double *A, *B, *RefOut;
//__device__ int bNumx,bNumy;
#define bNumx (128)
#define bNumy (128)
#define bSizex (N/bNumx)
#define bSizey (N/bNumy)
//__device__ int bSizex=N/bNumx;
//__device__ int bSizey=N/bNumy;


__global__ void p1_kernel(double*, double*, double*);

int main(){

    cudaMallocHost((void**)&A, N*N*sizeof(double));
    cudaMallocHost((void**)&B, N*N*sizeof(double));
    cudaMallocHost((void**)&RefOut, N*N*sizeof(double));

    double clkbegin, clkend;
    double t;

    int i,j,k;
    printf("Matrix Size = %d\n",N);
    printf("dimGrid(%d,%d), dimBlock(%d,%d)\n",bNumx,bNumy,bSizex,bSizey);

    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
            A[i*N+j] = 1.0*(i+0.25*j)/(N+1);
            B[i*N+j] = 1.0*(i-0.5*j)/(N+1);
        }


    clkbegin = rtclock();
    //
    // Time the reference baseline version of code 
    // Its output is to be compared to the test version to verify correctness
    //
    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
            for (k=i; k<N; k++)
                B[i*N+j] += A[i*N+k]*B[k*N+j];
    //
    // End of reference code
    //
    clkend = rtclock();
    t = clkend-clkbegin;
    if (B[(N/2)*N+(N/2)]*B[(N/2)*N+(N/2)] < -1000.0) 
        printf("To foil dead-code elimination by compiler: should never get here\n");
    printf("Base Sequential Symm-MatMult: %.1f GFLOPS; Time = %.3f sec; \n",
            1.0*N*N*(N+1)/t/1.0e9,t);

    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            RefOut[i*N+j] = B[i*N+j];

    // Initialize arrays

    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
            A[i*N+j] = 1.0*(i+0.25*j)/(N+1);
            B[i*N+j] = 1.0*(i-0.5*j)/(N+1);
        }

    double *d_A, *d_B, *d_C;
    cudaMalloc((void**)&d_A, N*N*sizeof(double));
    cudaMalloc((void**)&d_B, N*N*sizeof(double));
    cudaMalloc((void**)&d_C, N*N*sizeof(double));

    //
    // Test version of code; To be modified by you to improve performance
    //
    clkbegin = rtclock();

    // Send data to the GPU
    cudaMemcpy(d_A, A, N*N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, N*N*sizeof(double), cudaMemcpyHostToDevice);
    
    // Set grid and block sizes
    dim3 dimGrid(bNumx,bNumy);
    dim3 dimBlock(bSizex,bSizey);

    // Call the kernel (located below, you must modify it too)
    p1_kernel<<< dimGrid, dimBlock >>>(d_A, d_B, d_C);

    // Retrieve results
    cudaMemcpy(B, d_C, N*N*sizeof(double), cudaMemcpyDeviceToHost);

    clkend = rtclock();

    t = clkend-clkbegin;
    if (B[(N/2)*N+(N/2)]*B[(N/2)*N+(N/2)] < -1000.0) 
        printf("To foil dead-code elimination by compiler: should never get here\n");
    printf ("%.1f GFLOPS; Time = %.3f sec; \n",
            1.0*N*N*(N+1)/t/1.0e9,t);

    //
    // Verify correctness by comparing result with reference version's
    //
    compare();
}

//
// Test version of the kernel; To be modified by you to improve performance
// 
__global__ void p1_kernel(double *A, double *B, double *C) {
  //unsigned int i, j, k;
    //    for (i=0; i<N; i++)
    //  for (j=0; j<N; j++)
    //      for (k=i; k<N; k++)
    //          B[i*N+j] += A[i*N+k]*B[k*N+j];
    int bx=blockIdx.x; int by=blockIdx.y;
    int tx=threadIdx.x; int ty=threadIdx.y;
    int bSizeX= bSizex;
    int bSizeY= bSizey;
    
    int i = bSizeY*by+ty;
    int j = bSizeX*bx+tx;

    double Bsub=B[i*N+j];
    for (int k = i; k<N; k++)
      Bsub += A[i*N+k]*B[k*N+j];
    
    //_syncthreads();

    C[i*N+j]=Bsub;
    //printf("bx:%d by:%d tx:%d ty:%d i:%d j:%d\n", bx,by,tx,ty,i,j);
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
            this_diff = RefOut[i*N+j]-B[i*N+j];
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

