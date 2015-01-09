#include <omp.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#define N (128)
#define threshold (0.000000001)
double rtclock();
void compare();
double A[N][N][N], B[N][N][N], C[N][N], RefOut[N][N];

int main(){
    double clkbegin, clkend;
    double t;
    int nt, maxthr;

    int i,j,k,l;
    printf("Matrix Size = %d\n",N);

    for(i=0;i<N;i++)
        for(j=0;j<N;j++) {
            C[i][j] = 1.0*(i+0.25*j)/(N+1);
            for(k=0;k<N;k++)
            {
                A[i][j][k] = 1.0*(i+0.25*j-0.125*k)/(N+1);
                B[i][j][k] = 1.0*(i-0.5*j+0.125*k)/(N+1);
            }
        }


    clkbegin = rtclock();
    //
    // Time the reference baseline version of code 
    // Its output is to be compared to the test version to verify correctness
    //
    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
            for (k=0; k<N; k++)
                for (l=k+1; l<N; l++)
                    C[k][l] += 0.5*(A[l][i][j]*B[k][i][j]+A[k][i][j]*B[l][i][j]);
    //
    // End of reference code
    //
    clkend = rtclock();
    t = clkend-clkbegin;
    if (C[N/2][N/2]*C[N/2][N/2] < -1000.0) 
        printf("To foil dead-code elimination by compiler: should never get here\n");
    printf("Base Sequential: %.1f GFLOPS; Time = %.3f sec; \n",
            5.0*N*N*N*(N-1)/2.0/t/1.0e9,t);


    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            RefOut[i][j] = C[i][j];

    maxthr = omp_get_max_threads();
    printf("Maximum threads allowed by system is: %d\n",maxthr);

    // Loop to run test on diffreent number of threads

    for (nt=1;nt<=maxthr;nt++)
    {
        omp_set_num_threads(nt);

        // Initialize arrays

        for(i=0;i<N;i++)
            for(j=0;j<N;j++) {
                C[i][j] = 1.0*(i+0.25*j)/(N+1);
                for(k=0;k<N;k++)
                {
                    A[i][j][k] = 1.0*(i+0.25*j-0.125*k)/(N+1);
                    B[i][j][k] = 1.0*(i-0.5*j+0.125*k)/(N+1);
                }
            }

        printf("Requesting thrds=%d\n", nt);

        clkbegin = rtclock();
	int ii,jj;
#pragma omp parallel private(i,j,k,l,ii,jj)
        // Template version is intentionally made sequential
        // Make suitable changes to create parallel version
        {
            if (omp_get_thread_num()==0 && omp_get_num_threads() != nt)
                printf("Warning: Actual #threads %d differs from requested number %d\n",omp_get_num_threads(),nt);

            //
            // Test version of code; initially just contains a copy of base code
            // To be modified by you to improve performance
            //
	    //#pragma omp master
#pragma omp for
	    //	    for(ii=0;ii<N;ii+=N/2)
	    // for(jj=0;jj<N;jj+=N/2)
		for(k=0;k<N;k++)
		  for(l=k+1;l<N;l++)
		    for(ii=0;ii<N;ii+=N/2)
		      for(jj=0;jj<N;jj+=N/2)
			for (i=ii; i<ii+N/2; i++)
			  for (j=jj; j<jj+N/2; j++)
			    C[k][l] += 0.5*(A[l][i][j]*B[k][i][j]+A[k][i][j]*B[l][i][j]);
		
        } // End of parallel region
        clkend = rtclock();
        t = clkend-clkbegin;
        if (C[N/2][N/2]*C[N/2][N/2] < -1000.0) 
            printf("To foil dead-code elimination by compiler: should never get here\n");
        printf ("%.1f GFLOPS with %d threads; Time = %.3f sec; \n",
                5.0*N*N*N*(N-1)/2.0/t/1.0e9,nt,t);

        //
        // Verify correctness by comparing result with reference version's
        //
        compare();
    } // End of loop over different number of threads
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

