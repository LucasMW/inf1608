#include<stdio.h>
#include<stdlib.h>
#include<math.h>
static void MCopy(double * org, double* dest,int n)
{
	int i,n2;
	n2=n*n;
	for(i=0;i<n2;i++)
		dest[i]=org[i];
}
static void VCopy(double * org, double* dest,int n)
{
	int i;
	for(i=0;i<n;i++)
		dest[i]=org[i];
}
int MIndex(int i, int j,int n)
{
	return i*n+j;
}
void MPrint(double* A,int n)
{
	int i,j;
	for(i=0;i<n;i++)
	{	
		printf("|  ");
		for(j=0;j<n;j++)
		{
			printf("%.5lf  ",A[MIndex(i,j,n)]);
		}
		printf("|\n");
	}
	printf("\n");
}
void VPrint(double* V, int n)
{
	int i;
	printf("{ ");
	for(i=0;i<n;i++)
		printf("%lf ",V[i]);
	printf("}\n");
}
void Cholesky (int n, double* A)
{
	int k,j,i;
	double* ut =(double*)malloc(sizeof(double)*(n-1));
	//VCopy(A+1,ut,n-1); //copy first line
	for(k=0;k<n;k++)
	{
		A[MIndex(k,k,n)]=sqrt(A[MIndex(k,k,n)]);
		for(j=0;j<n-1;j++)
		{
			ut[j]=A[MIndex(j,j+1,n)]/(1.0*A[MIndex(j,j,n)]);
		}
		for(i=k+1;i<n;i++)
		{
			for(j=k+1;j<n;j++)
			{
				//A[MIndex(i,j,n)]=
			}
		}
	}
	
}

double Mult (int n,double* M, double*v, double* r)
{
	int i,j;
	double soma;
	for(i=0;i<n;i++)
	{
		soma=0;
		for(j=0;j<n;j++)
		{
			soma+=v[j]*M[MIndex(i,j,n)];
		}
		r[i]=soma;
	}
}
double Dot(int n,double*v,double* w)
{
	double prod=0;
	int i;
	for(i=0;i<n;i++)
	{
		prod+=v[i]*w[i];
	}
	return prod;
}
void ConjugateGradient(int n,double* A,double*b,double* x)
{
	
	int k;
	double* d = (double*) malloc(sizeof(double)*n);
	d[0];
	for(k=0;k<n-1;k++)
	{
	}
}
int main (void)
{
	return 0;
}
