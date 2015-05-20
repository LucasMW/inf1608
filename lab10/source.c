#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#define DEBUG //uncomment this line if you want to see the computing proccess
double Euler(double t0, double t1, double h, double y0, double (*f) (double t, double y))
{
	double y;
	double t;
	//t1=t1+h/2.0; //interval security
	for(t=t0,y=y0;t<t1;t+=h)
	{
		y=y+h*f(t,y);
	}
	return y;
}
double MidPoint(double t0, double t1, double h, double y0, double (*f) (double t, double y))
{
	double s1;
	double y;
	double t;
	//t1=t1+h/2.0; //interval security
	for(t=t0,y=y0;t<t1;t+=h)
	{
		s1=h*f(t,y);
		y=y+h*f(t+h/2,y+s1/2);
	}
	return y;
	
}
double RungeKutta(double t0, double t1, double h, double y0, double (*f) (double t, double y))
{
	double s1,s2,s3,s4;
	double y;
	double t;
	//t1=t1+h/2.0; //interval security
	for(t=t0,y=y0;t<t1;t+=h)
	{
		s1=h*f(t,y);
		s2=h*f(t+h/2,y+s1/2.0);
		s3=h*f(t+h/2,y+s2/2.0);
		s4=h*f(t+h,y+s3);
		y = y + (s1+2*(s2+s3)+s4)/6.0;
	}
	return y;
	
}
double FunctionDerivative(double t, double y)
{
	return (y+t*t)*t;
}
double FunctionRight(double t)
{
	return exp(t*t/2.0)-t*t-2;
}

int main (void)
{
	double error;
	double r,y;
	double h;
	double t=2.4;
	double t0=0;
	double y0=-1.0;
	int i;
	y=FunctionRight(t);
	printf("metodo\t\t h\t\t y(2.4)\t\t erro\n");
	for(h=0.1;h>=0.001;h/=10.0)
	{
		printf("-------\n");
		
		r=Euler(t0,t,h,y0,FunctionDerivative);
		error=(r-y)/y;
		printf("Euler\t\t %lf\t %lf\t %lf\n",h,r,error);
		r=MidPoint(t0,t,2*h,y0,FunctionDerivative);
		error=(r-y)/y;
		printf("MidPoint\t %lf\t %lf\t %lf\n",2*h,r,error);
		r=RungeKutta(t0,t,4*h,y0,FunctionDerivative);
		error=(r-y)/y;
		printf("RungeKutta\t %lf\t %lf\t %lf\n",4*h,r,error);
	}

	return 0;
}
