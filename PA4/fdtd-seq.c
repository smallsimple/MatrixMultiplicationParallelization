#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#define nx 4000
#define ny 4000
#define tmax 100
#define coeff1 0.5
#define coeff2 0.7

double ex[nx][ny];
double ey[nx][ny];
double hz[nx][ny];
void check();
void dump();

int main(){
  double rtclock();
				
  int tt,i,j,nt;
  double clkbegin, clkend;
  double t, maxdiff;
  
// Initialize arrays
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      ex[i][j] = sin(i)*(1-sin(j));
    }
  }
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      ey[i][j] = cos(i)*(1-cos(j));
    }
  }
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      hz[i][j] = sin(i)*(1-cos(j));
    }
  }

  clkbegin = rtclock();
  for (tt=0; tt<tmax; tt++){
    for (j=0; j<ny; j++)
      ey[0][j] = tt;

     for (j=0; j<ny; j++) 
      for (i=1; i<nx; i++) 
       ey[i][j] = ey[i][j] - coeff1*(hz[i][j]-hz[i-1][j]);

     for (j=1; j<ny; j++) 
      for (i=0; i<nx; i++) 
       ex[i][j] = ex[i][j] - coeff1*(hz[i][j]-hz[i][j-1]);

     for (j=0; j<ny-1; j++) 
      for (i=0; i<nx-1; i++) 
       hz[i][j] =  hz[i][j] -
         coeff2*(ex[i][j+1]-ex[i][j]+ey[i+1][j]-ey[i][j]);
  }
  clkend = rtclock();
  t = clkend-clkbegin;
  printf ("Sequential GFLOPS: %.1f, Time: %.1f\n", 10.0*nx*ny*tmax/t/1e9,t);
  check();
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


