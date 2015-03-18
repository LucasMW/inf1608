#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
double bisection(double a, double b, int p, double (*f) (double x, void* data), void* data)
{
	double fA;
	double fB;
	double fC;
	double c;
	double tolerance;

	tolerance = 0.5 * pow(10.0,-1.0*p);
	while(1)
	{
		//error = abs(a-b)/POW(2,n+1) use this to calculate interaction number
		fA=f(a,data);
		fB=f(b,data);
		c=(a+b)/2.0;
		fC=f(c,data);
		printf("%.6lf c:%.6lf\n",fC,c);
		if(fC==0 || fabs(b-a) < tolerance)
		{ //found
			printf("found \n");
			return c;
		}
		if(fC*fA>0)
			a=c;
		else
			b=c;
	}
	


}
double sqrtNet(double x,void* data)
{
	double* z=(double*)data;
	return (x+*z/x)/2.0;
}
double square_root(double z, int p)
{
	return bisection(0.0,z,p,sqrtNet,&z);
}
void square_rootTest(int p)
{
	int i;
	double r,x,res;
	res=0.5 * pow(10.0,-1.0*p);
	srand(time(NULL)); //changes random seed
	for(i=0;i<1000;i++)
	{
	r=(double)rand()/RAND_MAX;

	x=-0.1+ r*1000.0/2.0;
	if(!fabs(square_root(x,p)-sqrt(x))<(res))
	{
	printf("Error!\n");
	break;
	}
	}
	printf("All right\n");
	}
double Function(double x,void* data)
{
	return x*x*x - x*sin(2*x) + 0.2 ;
}
int main (void)
{
	double r1,r2,r3;
	printf("part1\n");
	 r1=bisection(-0.5,0.25,6,Function,NULL);
	 r2=bisection(0.25,0.5,6,Function,NULL);
	 r3=bisection(0.75,1.0,6,Function,NULL);
	printf("r1: %16g\nr2: %16g\nr3: %16g\n",r1,r2,r3);
	printf("part2\n");
	square_rootTest(6);
	 return 0;
}

