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

void Mult (int n,double* M, double*v, double* r)
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
int VetorIsZero(double* v ,int tam)
{
	int i;
	VPrint(v,tam);
	if(tam==0)
		return v[0];
	for(i=0;i<tam;i++)
	{
		if(v[i]!=0)
			return 0;
	}
	return 1;
}

void ConjugateGradient(int n,double* A,double*b,double* x)
{
	
	int k,i,j;
	double* d = (double*) malloc(sizeof(double)*n);
	double* r = (double*) malloc(sizeof(double)*n);
	double alpha;
	double beta;
	double* adk = (double*) malloc(sizeof(double)*n);
	Mult(n,A,x,r) ;
	for(i=0;i<n;i++)
	{
		r[i]= b[i] -r[i]; // r= b -Ax
	}
	VPrint(r,n);
	d[0]=r[0];
	for(k=0;k<n-1;k++)
	{
		printf("do\n");
		if(VetorIsZero(r,n))
		{	
				printf("is\n");
			break;
		}
		Mult(k,A,d,adk);
		alpha= Dot(k,r,r)/ (Dot(k,d,adk));
		printf("x: ");

		x[k+1]= x[k] + alpha * d[k];
		VPrint(x,k+1);
		beta= Dot(k,r+1,r+1) / Dot(k,r,r);
		d[k+1]= r[k+1] + beta * d[k];

	}
}
int main (void)
{
	double M1[] = 
	{	1.0,	-1.0,	0.0,
		-1.0,	2.0,	1.0,
		0.0,	1.0,	2.0
	};
	double M2[] = 
	{	1.0,	-1.0,	0.0,
		-1.0,	2.0,	1.0,
		0.0,	1.0,	5.0
	};
	double b1[] = {0,2,3};
	double b2[]= {3,-3,4};
	double x1[] ={0,0,0};
	double x2[]={1,1,1};
	ConjugateGradient(3,M1,b1,x1);
	ConjugateGradient(3,M2,b2,x2);
	VPrint(x1,3);
	VPrint(x2,3);
	return 0;
}
