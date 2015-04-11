#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif
//#define DEBUG //uncomment this line if you want to see the computing proccess
static long long int fact(int n)
{
	int i;
	long long int r=n;
	for(i=1;i<n;i++)
	{
		r*=(n-i);
	}
	return r;
}
double myCos(double x) //does not exist
{
	if(x<=M_PI/2)
	{
		 
	}
	else if(x>M_PI/2)
	{
		if(x>M_PI)
		{
			if(x>1.5*M_PI)
			{
				if(x>2*M_PI)
				{
					cos(x-fmod(x,2*M_PI));
				}
				else // ]3Pi/2, 2Pi]
				{
					return cos(2*M_PI-x);
				}
			}
			else // ]Pi. 3Pi/2]
			{
				return -myCos(x-M_PI);
			}
		}
		else // ]Pi/2,PI]
		{
			 return -myCos(M_PI-x);
		}
	}
	return 0.0;
}

int Chebyshev(int p, double ** X, double ** Y)
{
	double error;
	double tolerance;
	int n;
	int i;
	int b;
	tolerance=pow(10.0,-p);
	for(n=2;n<1000;n++)
	{
		error=(M_PI/2.0)/(pow(2.0,n-1)*fact(n));
		if(error<tolerance)
			break;
	}
	*X = (double*)malloc(sizeof(double)*n);
	*Y=  (double*)malloc(sizeof(double)*n);
	for(i=0,b=1;i<n;i++,b+=2)
	{
		
		(*X)[i]=(M_PI/2)/2.0*cos(b*M_PI)/(2*n)+(M_PI/2.0)/2.0;
		(*Y)[i]=cos((*X)[i]);
	}
	return n;
}
double Lagrange (int n,double* X, double* Y, double x)
{
	return 0.0;
}
double MyCos1(int n, double* X, double* Y,double x)
{
	return 0.0;
}
double* NewtonCoefficients(int n, double* X, double* Y)
{
}
double Newton(int n, double* X, double* Y,double* b,double x)
{
}
double MyCos2(int n, double* X, double* Y, double* b, double x)
{
	return 0.0;
}
int main (void)
{
	double *amostrasX, *amostrasY;
	int n;
	n=Chebyshev(6,&amostrasX,&amostrasY);
	printf("n: %d \n",n);
	printf("fact %d\n",(int)fact(5));
	return 0;
}
