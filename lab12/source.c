#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#define DEBUG //uncomment this line with want to see the computing proccess

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
int checkEDD(double* A, int n)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		double soma=0;
		for(j=0;j<n;j++)
		{
			if(i!=j)
				soma+=A[MIndex(i,j,n)];
		}
		if(fabs(A[MIndex(i,i,n)])<fabs(soma)) //Not Dominant
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

	if(!checkEDD(A,n))
	{
		printf("Not Diagonal Dominant\nExiting\n");
		return;
	}
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
					#ifdef DEBUG
					//MPrint(A,n);
					printf("i %d j %d  k %d s %lf\n",i,j,k,s);
					#endif
				}
			}
			x[i]=(b[i]-s)/A[MIndex(i,i,n)];

		}
		#ifdef DEBUG
		//sleep(1);
		 VPrint(x,n);
		#endif
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
	if(!checkEDD(A,n))
	{
		printf("Not Diagonal Dominant\nExiting\n");
		return;
	}
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
	if(!checkEDD(A,n))
	{
		printf("Not Diagonal Dominant\nExiting\n");
		return;
	}
	for(k=0;k<niter;k++)
	{
		for(i=0;i<n;i++)
		{
			s=0;
			for(j=0;j<n;j++)
			{
				if(j!=i)
					s+=A[MIndex(i,j,n)]*x[j];
				#ifdef DEBUG
					printf("i %d j %d  k %d s %lf\n",i,j,k,s);
					#endif
			}
			x[i]+=w*((b[i]-s)/A[MIndex(i,i,n)]-x[i]);
		}
		#ifdef DEBUG
		sleep(1);
		 VPrint(x,n);
		#endif
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
	double xr1[6];
	double M2[]=
	{
		1.0,	2.0,	-1.0,
		2.0,	1.0,	-2.0,
		-3.0,	1.0,	1.0
	};
	double b2[]= {3.0,3.0,-6.0}; 
	double n2=3;
	double x2[]= {0.0,0.0,0.0};
	double xr2[3];
	int steps= 5;
	printf("Matrix 1\n");
	MPrint(M1,n1);
	printf("Jacobi \n");
	Jacobi(n1,M1,b1,x1,steps);
	VPrint(x1,n1);
	VCopy(x1,xr1,n1);
	printf("GaussSeidel \n");
	GaussSeidel(n1,M1,b1,x1,steps);
	VPrint(x1,n1);
	VCopy(x1,xr1,n1);
	printf("SOR \n");
	SOR(n1,M1,b1,x1,steps,1.2);
	VPrint(x1,n1);
	VCopy(x1,xr1,n1);

	printf("Matrix 2\n");
	MPrint(M2,n2);
	printf("Jacobi\n");
	Jacobi(n2,M2,b2,x2,steps);
	VPrint(x2,n2);
	VCopy(x2,xr2,n2);
	printf("GaussSeidel \n");
	GaussSeidel(n2,M2,b2,x2,steps);
	VPrint(x2,n2);
	VCopy(x2,xr2,n2);
	printf("SOR \n");
	SOR(n2,M2,b2,x2,steps,1.2);
	VPrint(x2,n2);
	VCopy(x2,xr2,n2);

	return 0;
}
