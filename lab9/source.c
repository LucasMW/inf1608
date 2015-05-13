#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.1415926535
#endif
double Simpson (double a, double b, double (*f) (double x))
{
	double mid=(a+b)/2.0;
	double h = b-a;
	double result;

	result=f(a)+4*f(mid)+f(b);
	result*=h/6.0;
	return result;
}

double DoubleSimpson(double a, double b, double (*f) (double x))
{
	double mid=(a+b)/2.0;
	double h = b-a;
	double result;
	double error;
	result=Simpson(a,mid,f)+Simpson(mid,b,f);
	error=15.0*(Simpson(a,b,f)-Simpson(a,mid,f)-Simpson(mid,b,f));
	return result;
}
double Function1(double x)
{
	return x/sqrt(x*x*x*x+1);
}
double Function2(double x)
{
	return x*x*sin(x);
}
double AdaptativeSimpson(double a, double b, double (*f) (double x), double tol)
{

}

int main (void)
{
	printf("Simpson Tests\n");
	printf("int(f1,%lf,%lf)=%lf\n",0.0,1.0,Simpson(0.0,1.0,Function1));
	printf("int(f2,%lf,%lf)=%lf\n",0.0,M_PI,Simpson(0.0,M_PI,Function2));

	return 0;
}
