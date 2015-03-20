#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
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
		//printf("%.6lf c:%.6lf\n",fC,c);
		if(fC==0 || fabs(b-a) < tolerance)
		{ //found
			//printf("found \n");
			return c;
		}
		if(fC*fA>0)
			a=c;
		else
			b=c;
	}
	


}
double sqrtNet(double x,double z)
{
	
	return (x + z/x)/2.0;
}
double square_root(double z, int p)
{
	double x, gx;
	double tolerance;
	tolerance = 0.5 * pow(10.0,-1.0*p);
	x=z/2.0;
	gx=sqrtNet(x,z);
	while(fabs(gx-x)>tolerance)
	{	
		x=gx;
		gx=sqrtNet(x,z);
	}
	return x;
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

		x=0.0+ r*1000.0/2.0;
		if(!(fabs(square_root(x,p)-sqrt(x))<res))
		{
			printf("Error!\n");
			return;
		}
	}
	printf("All right\n");
}
double Function(double x,void* data)
{
	return x*x*x - x*sin(2*x) + 0.2 ;
}
double Heq(double h, void* data)
{
	double *coef;
	double h2;
	coef=(double*)data;
	h2=h*h;

	return -1.0*M_PI*h2*h + 3*M_PI*h2*coef[0] -3*coef[1];


}
double height(double r, double v)
{

	double vEsfera=M_PI*(4.0/3.0)*r*r*r;
	double coeficientes[2];
	if(v>vEsfera)
	{

		return -1;
	}
	coeficientes[0]=r;
	coeficientes[1]=v;
	return bisection(0.0,2*r,6,Heq,coeficientes);

}
int main (void)
{
	double r1,r2,r3;
	double h,r,v;
	printf("part1\n");
	 r1=bisection(-0.5,0.25,6,Function,NULL);
	 r2=bisection(0.25,0.5,6,Function,NULL);
	 r3=bisection(0.75,1.0,6,Function,NULL);
	printf("r1: %16g\nr2: %16g\nr3: %16g\n",r1,r2,r3);
	printf("part2\n");
	//printf(" %lf VS %lf",square_root(0.5,6),sqrt(0.5));
	square_rootTest(6);
	printf("part3\n");
	r=10.0;
	v=4000.0;
	h=height(10.0,4000.0);
	if(h<0)
	{
		printf("Erro: Altura nao pode ser calculada!\n");
	}
	printf("h:%lf r:%lf v:%lf \n", h,10.0,4000.0);

	 return 0;
}

