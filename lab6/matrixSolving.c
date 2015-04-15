
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


static void MCopy(double * org, double* dest,int n)
{
	int i,n2;
	n2=n*n;
	for(i=0;i<n2;i++)
		dest[i]=org[i];
	
	
}
static double MElement(double* A,int i, int j,int n)
{
	return *(A+(i*n+j));
}
static double VPrint(double* V, int n)
{
	int i;
	printf("{ ");
	for(i=0;i<n;i++)
		printf("%lf ",V[i]);
	printf("}\n");
}

int MIndex(int i, int j,int n)
{
	return i*n+j;
}
double MPrint(double* A,int n)
{
	int i,j;

	for(i=0;i<n;i++)
	{	
		printf("|  ");
		for(j=0;j<n;j++)
		{
			printf("%lf  ",A[MIndex(i,j,n)]);
		}
		printf("|\n");
		
	}
	printf("\n");
}
void NaiveGaussElimination(int n, double* A, double* b)
{
	//Ax=b : A is n*n Matrix
	int i,j,k;
	double f;
	for(j=0; j<n-1;j++)
	{
		for(i=j+1; i<n; i++)
		{
			f=MElement(A,i,j,n) / MElement(A,j,j,n);
			#ifdef DEBUG
			printf("f: %lf\n",f);
			#endif
			for(k=j; k<n; k++)
			{
				A[MIndex(i,k,n)]-=f*A[MIndex(j,k,n)];
				#ifdef DEBUG
				MPrint(A,n);
				#endif
			}
			b[i]-=f*b[j];
		}
		
	}
}
void BackSubstitution(int n, double* A, double* b, double* x)
{
	int i,j;
	double s;
	for(i=n-1;i>=0;i--)
	{
		s=0;
		for(j=i+1;j<n;j++)
		{
			s+=A[MIndex(i,j,n)]*x[j];
			#ifdef DEBUG
			printf("s:%lf j:%d i:%d\n",s,j,i);
			#endif
		}
		x[i]=(b[i]-s)/A[MIndex(i,i,n)];
		#ifdef DEBUG
		VPrint(x,n);
		#endif
	}
}
