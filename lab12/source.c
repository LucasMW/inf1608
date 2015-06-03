#include<stdio.h>
#include<stdlib.h>

#define DEBUG //uncomment this line with want to see the computing proccess

static void MCopy(double * org, double* dest,int n)
{
	int i,n2;
	n2=n*n;
	for(i=0;i<n2;i++)
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
int checkEDD(double* A, int n)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		double soma=0;
		for(j=0;j<n;i++)
		{
			if(i!=j)
				soma+=A[MIndex(i,j,n)];
		}
		if(A[MIndex(i,j,n)]<soma)
		{
			return 0;
		}
			
	}
	return 1;
}
void Jacobi (int n , double* A, double* b, double* x, int niter)
{
	double * L;
	double * D;
	double * U;
	int i,j;
	int k; //count iterations
	double s;
	for(k=0;k<niter;k++)
	{
	for(i=0;i<n;i++)
	{
		s=0;
		for(j=0;j<n;j++)
		{
			if(i!=j)
			{
				s+=A[MIndex(i,j,n)]*x[j];
			}
		}
		x[i]=(b[i]-s)/A[MIndex(i,j,n)];
	}
	}

	/*for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			if(i<j) //superior
				U[MIndex(i,j,n)]=A[MIndex(i,j,n)];
			else if(i>j) //inferiro
				L[MIndex(i,j,n)]=A[MIndex(i,j,n)];*/

}
void GaussSeidel (int n, double* A, double* b, double* x, int niter)
{
	int i,j,k;
	double soma;
	for(k=0; k<niter; k++)
	{
		for(i=0;i<n;i++)
		{
			soma=0;
			for(j=0;j<n;j++)
			{
				if(j!=i)
					soma+=A[MIndex(i,j,n)]*x[j];
			}
			x[i]=(b[i]-soma)/A[MIndex(i,i,n)];

		}

	}
}

void SOR(int n, double* A, double * b, double* x, int niter, double w)
{
	int i,j,k;
	double s;
	for(k=0;k<niter;k++)
	{
		for(i=0;i<niter;i++)
		{
			s=0;
			for(j=0;j<n;j++)
			{
				if(j!=i)
					s+=A[MIndex(i,j,n)]*x[j];
			}
			x[i]+=w*((b[i]-s)/A[MIndex(i,i,n)]-x[i]);
		}
	}
}
int main (void)
{
	double M1[] =
	{
		3.0,	-1.0,	0.0,	0.0,	0.0,	0.5,
		-1.0,	3.0,	-1.0,	0.0,	0.5,	0,
		0.0,	-1.0,	3.0,	-1.0,	0.0,	0.0,
		0.0,	0.0,	-1.0,	3.0,	-1.0,	0.0,
		0.0,	0.5,	0.0,	-1.0,	3.0,	-1.0,
		0.5,	0.0,	0.0,	0.0,	-1.0,	3.0
	};
	double b1[]= {2.5,1.5,1.0,1.0,1.5,2.5};
	double n1=6;
	double x1[] ={0.0,0.0,0.0,0.0,0.0,0.0};
	double M2[]=
	{
		1.0,	2.0,	-1.0,
		2.0,	1.0,	-2.0,
		-3.0,	1.0,	1.0
	};
	double b2[]= {3.0,3.0,-6.0}; 
	double n2=6;
	double x2[]= {0.0,0.0,0.0};
	Jacobi(n1,M1,b1,x1,5);
	VPrint(x1,n1);
	GaussSeidel(n1,M1,b1,x1,5);
	VPrint(x1,n1);
	SOR(n1,M1,b1,x1,5,1.2);
	VPrint(x1,n1);

	return 0;
}
