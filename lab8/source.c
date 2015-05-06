#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double functionDerivative(double x)
{
	return cos(x)+sin(x);
}
double function(double x)
{
	return sin(x)-cos(x);
}
double Error(double h)
{
	double erroMaq=9.2-9.0-0.2;

	return h*h/6.0 *1.415 + erroMaq/h;
}
double Derivative(double x, double h,double(*f)(double x))
{
	return (f(x+h)-f(x-h))/(2*h);
}
int main (void)
{
	int i;
	double h;
	double r;

	for(i=1;i<=12;i++)
	{
		h=pow(10.0,(double)-i);
		r=Derivative(0.0,h,function);
		printf("%16g %16g %16g\n",h,functionDerivative(0)-r);
	}
	return 0;
}
