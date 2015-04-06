#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#define DEBUG //uncomment this line to see the computing proccess
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

static int MIndex(int i, int j,int n)
{
	return i*n+j;
}
static double MPrint(double* A,int n)
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
void GaussElimination(int n, double* A, double* b)
{ 
	int i,j,k;
	int index;
	char flag;
	double pvt,t;
	double f;
	for(j=0; j<n-1;j++)
	{
		for(i=j+1; i<n; i++)
		{
			//Finding Pivot
			pvt=A[MIndex(i,i,n)];
			index=i;
			flag=0;
			for(k=i+1;k<n;k++)
			{
				if(pvt<A[MIndex(k,i,n)])
				{	printf("pvt: %lf\n",pvt);
					index=k;
					flag=1;
					pvt=fabs(A[MIndex(k,i,n)]);
				}

			}
			//printf("%16g",pvt-0);
			//	return;// Matrix is singular
			
			//Swap Rows
		if(flag)
		{
			for(k=0;k<n;k++)
			{
				t=A[MIndex(index,k,n)];
				A[MIndex(index,k,n)]=A[MIndex(i,k,n)];
				A[MIndex(i,k,n)]=t;
			}
			#ifdef DEBUG
			printf("Swapped\n");
			MPrint(A,n);
			#endif
		}
			
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
	double xV2[6];
	double M3[2*2]={pow(10,-17),1.0,1.0,2.0};
	double b3[2]={1.0,4.0};
	double xV3[2];
	double M1c[3*3],M2c[6*6],M3c[2*2];
	double b1c[3]={3,3,-6},b2c[6]={2.5,1.5,1,1,1.5,2.5},b3c[2]={1.0,4.0};;
	int i;
	
	
	printf("part1\n");
	MCopy(M1,M1c,3);
	MCopy(M2,M2c,6);
	MCopy(M3,M3c,2);
	
	NaiveGaussElimination(3,M1,b1);
	BackSubstitution(3,M1,b1,xV1);
	for(i=0;i<3;i++)
	{
		printf("x%d:%16g\n",i,xV1[i]);
	}
	
	NaiveGaussElimination(6,M2,b2);
	printf("subs\n");
	MPrint(M2,6);
	BackSubstitution(6,M2,b2,xV2);
	for(i=0;i<6;i++)
	{
		printf("x%d:%16g\n",i,xV2[i]);
	}
	
	NaiveGaussElimination(2,M3,b3);
	printf("subs\n");
	BackSubstitution(2,M3,b3,xV3);
	for(i=0;i<2;i++)
	{
		printf("x%d:%16g\n",i,xV3[i]);
	}
	printf("part2\n");
	
	
	GaussElimination(3,M1c,b1c);
	printf("subs\n");
	BackSubstitution(3,M1c,b1c,xV1);
	for(i=0;i<3;i++)
	{
		printf("x%d:%16g\n",i,xV1[i]);
	}
	
	GaussElimination(6,M2c,b2c);
	printf("subs\n");
	BackSubstitution(6,M2c,b2c,xV2);
	for(i=0;i<6;i++)
	{
		printf("x%d:%16g\n",i,xV2[i]);
	}
	printf("subs\n");
	GaussElimination(2,M3c,b3c);
	BackSubstitution(2,M3c,b3c,xV3);
	for(i=0;i<2;i++)
	{
		printf("x%d:%16g\n",i,xV3[i]);
	}
	return 0;
}
