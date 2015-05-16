#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
double Simpson (double a, double b, double (*f) (double x))
{
	double mid=(a+b)/2.0;
	double h = b-a;
	double result;
	double error;
	result=f(a)+4*f(mid)+f(b);
	result*=h/6.0;
	error= h/2.0;
	error*=error;
	error*=error;
	error/=-90.0;
	printf("error is %lf\n",error)
;	return result;
}

double DoubleSimpson(double a, double b, double (*f) (double x))
{
	double mid=(a+b)/2.0;
	double h = b-a;
	double result;
	double error;
	result=Simpson(a,mid,f)+Simpson(mid,b,f);
	error=15.0*(Simpson(a,b,f)-Simpson(a,mid,f)-Simpson(mid,b,f));
	printf("error is %lf\n",error);
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
	double S;
	double E;
	double D;
	double h;

	h = b-a;

	S=DoubleSimpson(a,b,f);
	D=DoubleSimpson(a,a+h/2,f)+DoubleSimpson(a+h/2,b,f);
	E = (S-D)/15.0;
	sleep(1);
	printf("a:%lf b:%lf tol:%lf\n",a,b,tol);	
	if(fabs(E)>tol)
	{
		printf("Entered Tol S %lf E %lf D %lf tol %lf E-tol %lf\n",S,E,D,tol,E-tol);
		return AdaptativeSimpson(a,a+h/2.0,f,tol/2.0)+AdaptativeSimpson(a+h/2.0,b,f,tol/2.0);
	}
	else
	{
		printf("Passed S %lf E %lf D %lf tol %lf\n",S,E,D,tol);	
		return D;
	}

}

int main (void)
{
	printf("Simpson Tests\n");
	printf("int(f1,%lf,%lf)=%lf\n",0.0,1.0,Simpson(0.0,1.0,Function1));
	printf("int(f2,%lf,%lf)=%lf\n",0.0,M_PI,Simpson(0.0,M_PI,Function2));
	printf("Double Simpson Tests\n");
	printf("int(f1,%lf,%lf)=%lf\n",0.0,1.0,DoubleSimpson(0.0,1.0,Function1));
	printf("int(f2,%lf,%lf)=%lf\n",0.0,M_PI,DoubleSimpson(0.0,M_PI,Function2));
	printf("Adaptative Simpson Tests\n");
	printf("int(f1,%lf,%lf)=%lf\n",0.0,1.0,AdaptativeSimpson(0.0,1.0,Function1,pow(10.0,-5.0)));
	printf("int(f2,%lf,%lf)=%lf\n",0.0,M_PI,AdaptativeSimpson(0.0,M_PI,Function2,pow(10.0,-5.0)));

	return 0;
}
