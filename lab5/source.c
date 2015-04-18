#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif
//#define DEBUG //uncomment this line if you want to see the computing proccess

static void printVector(double* vect, int tam)
{
	int i;
	printf("< ");
	for(i=0;i<tam;i++)
		printf("%lf ",vect[i]);
	printf(">\n");
}


static long long int fact(unsigned int n)
{
	unsigned int i;
	long long int r=1;
	for(i=0;i<n;i++)
	{
		r*=(n-i);
	}
	return r;
}


int Chebyshev(int p, double ** X, double ** Y)
{
	double error;
	double tolerance;
	int n;
	int i;
	int b;
	double start, end; //the exist make code more legible and useful in the future
	double diffBy2,sumBy2,PiBy2n; //constant coeficients
	double value;
	tolerance=pow(10.0,-p);
	for(n=2;n<1000;n++)
	{
		error=(M_PI/2.0)/(pow(2.0,n-1)*fact(n));
		if(error<tolerance)
			break;
	}
	printf("n: %d with Error: %16g\n",n,error);
	start=0.0;
	end=M_PI/2.0;

	diffBy2=(end-start)/2.0;
	sumBy2=(start+end)/2.0;
	PiBy2n=M_PI/(n*2.0);
	*X=(double*)malloc(sizeof(double)*n);
	*Y=(double*)malloc(sizeof(double)*n);
	for(i=0,b=1;i<n;i++,b+=2)
	{
		value= diffBy2 * cos(b*PiBy2n)+ sumBy2;  //runs faster
		(*X)[i]=value;
		#ifdef DEBUG
		printf("X[%d]: %lf with b:%d ",i,value,b);
		printf("X[i] = %16g\n",(*X)[i]);
		#endif
		(*Y)[i]=cos((*X)[i]);
	}
	return n;
}
double Lagrange (const int n,double* X, double* Y, double x)
{
	double* L;
	int i,j;
	double result;

	L=(double*)malloc(sizeof(double)*n);
	for(i=0;i<n;i++)
	{
		L[i]=1;
		for(j=0;j<n;j++)
		{
			if(j!=i) // there's no (x-X[j])/(X[j]-X[j])
			{
				L[i]*=(x-X[j])/(X[i]-X[j]);
			}
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

				if(x>2.0*M_PI)
				{	
					#ifdef DEBUG
					printf("x:%lf e ]2PI, INF]\n",x);
					#endif
					return MyCos1(n,X,Y,fmod(x,2.0*M_PI));
				}
				else // ]3Pi/2, 2Pi]
				{
					#ifdef DEBUG
					printf("x:%lf e ]3Pi/2, 2Pi]\n",x);
					#endif
					return MyCos1(n,X,Y,2.0*M_PI-x);
				}
			}
			else // ]Pi, 3Pi/2]
			{
				#ifdef DEBUG
				printf("x:%lf e ]Pi, 3Pi/2]\n",x);
				#endif
				return -MyCos1(n,X,Y,x-M_PI);
			}
		}
		else // ]Pi/2,PI]
		{
			#ifdef DEBUG
			printf("x:%lf e ]Pi/2,PI]\n",x);
			#endif
			return -MyCos1(n,X,Y,M_PI-x);
		}
	}
	return 0.0;
	
}
void MyCos1Test(int p)
{
	double *amostrasX, *amostrasY;
	int n;
	int i;
	double r,x,res;
	double fTest,fControl,cmp;
	res=0.5 * pow(10.0,-1.0*p);
	srand(time(NULL)); //changes random seed
	n=Chebyshev(p,&amostrasX,&amostrasY);
	for(i=0;i<1000;i++)
	{
		r=(double)rand()/RAND_MAX;
		x=0.0+ r*20.0*M_PI;
		fTest=MyCos1(n,amostrasX,amostrasY,x);
		fControl=cos(x);
		cmp=fabs(fTest-fControl);
		if(!(cmp<res))
		{
			printf("Error in turn %d! \n%lf VS %lf \nerror: %lf  >= res %lf\nvalue: %lf\n",i+1,fTest,fControl,cmp,res,x);
			return;
		}
	}
	printf("All right\n");
	
}
// fdiv recursive impletation where i>j
double fdiv(int i, int j,double* X, double* Y)
{
	if(i==j)
	{
		return Y[i];
	}
	else
		return (fdiv(i,j+1,X,Y)-fdiv(i-1,j,X,Y))/(X[i]-X[j]);
}
double* NewtonCoefficients(int n, double* X, double* Y)
{
	double* b; // Newton coeficients
	int i,j;
	b=(double*)malloc(sizeof(double)*n); //allocating coefficients vector

	for(i=0;i<n;i++)
	{
		b[i]=fdiv(i,0,X,Y); //populating vector
	}
	return b;

}
double Newton(int n, double* X, double* Y,double* b,double x)
{
	double result;
	int i;
	result=b[n-1];
	for(i=n-2;i>=0;i--)
	{
		result *= (x-X[i]);
		result+=b[i];
	}
	return result;
}
double MyCos2(int n, double* X, double* Y, double* b, double x)
{
	if(x<=M_PI/2)
	{
		return Newton(n,X,Y,b,x); 
	}
	else if(x>M_PI/2)
	{
		if(x>M_PI)
		{
			if(x>1.5*M_PI)
			{

				if(x>2.0*M_PI)
				{
					#ifdef DEBUG
					printf("x:%lf e ]2PI, INF]\n",x);
					#endif
					return MyCos2(n,X,Y,b,fmod(x,2.0*M_PI));
				}
				else // ]3Pi/2, 2Pi]
				{
					#ifdef DEBUG
					printf("x:%lf e ]3Pi/2, 2Pi]\n",x);
					#endif
					return MyCos2(n,X,Y,b,2.0*M_PI-x);
				}
			}
			else // ]Pi, 3Pi/2]
			{
				#ifdef DEBUG
				printf("x:%lf e ]Pi, 3Pi/2]\n",x);
				#endif
				return -MyCos2(n,X,Y,b,x-M_PI);
			}
		}
		else // ]Pi/2,PI]
		{
			#ifdef DEBUG
			printf("x:%lf e ]Pi/2,PI]\n",x);
			#endif
			return -MyCos2(n,X,Y,b,M_PI-x);
		}
	}
	return 0.0;
}
void MyCos2Test(int p)
{
	double *amostrasX, *amostrasY;
	int n;
	int i;
	double* b;
	double r,x,res;
	double fTest,fControl,cmp;
	res=0.5 * pow(10.0,-1.0*p);
	srand(time(NULL)); //changes random seed
	n=Chebyshev(p,&amostrasX,&amostrasY);
	b=NewtonCoefficients(n,amostrasX,amostrasY);

	for(i=0;i<1000;i++)
	{
		r=(double)rand()/RAND_MAX;
		x=0.0+ r*20.0*M_PI;
		fTest=MyCos2(n,amostrasX,amostrasY,b,x);
		fControl=cos(x);
		cmp=fabs(fTest-fControl);
		if(!(cmp<res))
		{
			printf("Error in turn %d! \n%lf VS %lf \nerror: %lf  >= res %lf\nvalue: %lf\n",i+1,fTest,fControl,cmp,res,x);
			return;
		}
	}
	printf("All right\n");
	
}
int main (void)
{
	
	printf("Part1\n");
	MyCos1Test(10);
	printf("Part2\n");
	MyCos2Test(10);
	
	return 0;
}
