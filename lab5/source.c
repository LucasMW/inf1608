#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif
//#define DEBUG //uncomment this line if you want to see the computing proccess
static long long int fact(unsigned int n)
{
	int i;
	long long int r=1;
	for(i=0;i<n;i++)
	{
		r*=(n-i);
	}
	return r;
}
double myCos(double x) //does not exist
{
	
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
	double L[n];
	int i,j;
	double result;
	for(i=0;i<n;i++)
	{
		L[i]=1;
		for(j=0;j<n;j++)
		{
			if(j!=i) // there's no (x-X[j])/(X[j]-X[j])
			L[i]*=(x-X[j])/(X[i]-X[j]);
		}
	}	
	result=0.0;
	for(i=0;i<n;i++)
	{
		result+=(Y[i]*L[i]);
	}	

	return result;
}
double MyCos1(int n, double* X, double* Y,double x)
{
	 if(x<=M_PI/2)
	{
		return Lagrange(n,X,Y,x); 
	}
	else if(x>M_PI/2)
	{
		if(x>M_PI)
		{
			if(x>1.5*M_PI)
			{

				if(x>2*M_PI)
				{
					MyCos1(n,X,Y,x-fmod(x,2*M_PI));
				}
				else // ]3Pi/2, 2Pi]
				{
					return MyCos1(n,X,Y,2*M_PI-x);
				}
			}
			else // ]Pi. 3Pi/2]
			{
				return -MyCos1(n,X,Y,x-M_PI);
			}
		}
		else // ]Pi/2,PI]
		{
			 return -MyCos1(n,X,Y,M_PI-x);
		}
	}
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
	n=Chebyshev(10,&amostrasX,&amostrasY);
	
	printf("n: %d \n",n);
	printf("fact %d\n",(int)fact(4));
	printf("cos %16g\n",MyCos1(n,amostrasX,amostrasY,M_PI/4));
	return 0;
}
