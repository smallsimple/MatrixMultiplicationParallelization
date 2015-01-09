#include <stdio.h>
#include <stdint.h>
#include <math.h>
int main(int argc, char *argv[])
{
        long long int  n, size;
	int id,p;
	long long int low_value, high_value;
	n=atoll(argv[1]);
	p=atoll(argv[2]);
	printf("n=%d  p=%d\n",n,p);
	for (id = 0; id < p; id++)
	{
		low_value= 2 + id*(n-1)/p;
		//printf("low_value=%d\n",low_value);
		high_value=1 + (id+1)*(n-1)/p;
		//printf("high_value=%d\n",high_value);
		size = high_value - low_value +1;
		printf("id=%2lld low_value=%10lld high_value=%10lld size=%10lld\n",id,low_value,high_value,size);
		//printf("2+id*(n-1)/p=%d\n",2+id*(n-1)/p);
	}
}
