#include <stdio.h>
#include <stdlib.h>

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

static int MIndex(int i, int j,int n)
{
	return i*n+j;
}
static double MPrint(double* A,int n)
{
	int i,j;

	for(i=0;i<n;i++)
	{	
		printf("|\t");
		for(j=0;j<n;j++)
		{
			printf("%lf\t",A[MIndex(i,j,n)]);
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
	for(j=0; j<n;j++)
	{
		for(i=j+1; i<n; i++)
		{
			f=MElement(A,i,j,n) / MElement(A,j,j,n);
			printf("f: %lf\n",f);
			for(k=j; k<n; k++)
			{
				A[MIndex(i,k,n)]-=f*A[MIndex(j,k,n)];
				MPrint(A,n);
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
		for(j=i;j<n;j++)
		{
			s+=A[MIndex(i,j,n)]*x[j];
			printf("s:%lf j:%d i:%d\n",s,j,i);
		}
		x[i]=(b[i]-s)/A[MIndex(i,i,n)];
		VPrint(x,n);
	}
}
int main (void)
{
	double M1[3*3]={1,2,-1,2,1,-2,-3,1,1};
	double xV1[3]={0,0,0};
	double b1[3]={3,3,-6};
	double M2[6*6]=
	{	3,  -1,   0,  0,  0,   0.5,
	   	-1,  3,  -1,  0,  0.5, 0,
		0,  -1,   3, -1,  0,   0,
		0,   0,  -1,  3, -1,   0,
		0,   0.5, 0, -1,  3,  -1,
		0.5, 0,   0,  0, -1,   3
	};
	double b2[6]={2.5,1.5,1,1,1.5,2.5};
	double xV2[3];
	int i;
	
	printf("part1\n");
	NaiveGaussElimination(3,M1,b1);
	printf("subs\n");
	BackSubstitution(3,M1,b1,xV1);
	for(i=0;i<3;i++)
	{
		printf("%16g\n",xV1[i]);
	}
	

	return 0;
}
