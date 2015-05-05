#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "matrixSolving.h"

static int index(int i,int j, int n)
{
	return i*n+j;
}
double* cria_matriz(int m,int n)
{
	return (double*)malloc(sizeof(double)*m*n);
}
static void copia_matriz(int m,int n, double* org, double* dest)
{
	int i;
	int tam;
	tam=m*n;
	for(i=0;i<tam;i++)
		dest[i]=org[i];
}
static void print_matriz(int m,int n,double* A)
{
	int i,j;

	for(i=0;i<m;i++)
	{	
		printf("| ");
		for(j=0;j<n;j++)
			printf("%lf ",A[index(i,j,n)]);
		printf("|\n");
	}
}
static void print_vetor(double* v, int tam)
{
	int i;
	printf("{ ");
	for(i=0;i<tam;i++)
		printf("%lf ",v[i]);
	printf("}\n");
}
double* cria_transposta(int m,int n,double* A)
{
	int i;
	int j;
	double* T;
	T=cria_matriz(n,m);
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			T[index(j,i,m)]=A[index(i,j,n)];
	return T;
}
double* cria_mult(int m,int n,int l, double* A, double* B)
{				//A (mxn) * B (nxl)
	double * C;
	int i,j,k;
	double sum;
	C=cria_matriz(m,l);

	print_matriz(m,n,A);
	print_matriz(n,l,B);
	for(i=0;i<m;i++)
	{
		for(j=0;j<l;j++)
		{	sum=0;
			for(k=0;k<n;k++)
			{
				sum+=A[index(i,k,n)]*B[index(k,j,l)];
				printf("A(%d,%d) * B(%d,%d) \n%lf*%lf\nsum: %lf\n",i,k,k,j,A[index(i,k,n)],B[index(k,j,l)],sum);
			}
			printf("(%d,%d) %lf \n",i,j,sum);
			C[index(i,j,l)]=sum;
		}
	}
			return C;
}

double* cria_mmq(int m,int n,double* A,double* b)
{
	double* At;
	double* AtA;
	double* Atb;

	At=cria_transposta(m,n,A);
	AtA=cria_mult(n,m,n,At,A);
	Atb=cria_mult(n,m,1,At,b);

	return NULL;
}
double norma2(int n,double* v)
{
	int i;
	double soma;
	for(i=0,soma=0;i<n;i++)
		soma+=v[i]*v[i];
	return sqrt(soma);
}
int main (void)
{
	double *M1,*M2,*M3,*M4,*M5;
	double A[] = 
	{
		1.0,2.0,3.0,
		4.0,5.0,6.0,
		7.0,8.0,9.0
	};
	double B[] = 
	{
		1.0,2.0,
		4.0,5.0,
		7.0,8.0
	};
	double C[] = 
	{
		2.0,3.0,
		0.0,1.0,
		-1.0,4.0
	};
	double D[] = 
	{
		1.0,2.0,3.0
		-2.0,0.0,4.0
	};
	double E[] = 
	{
		1.0,2.0,
		3.0,4.0
	};
	double F[] = 
	{
		-1.0,-2.0,
		-3.0,-4.0
	};
	double I2[] = 
	{
		1.0,0.0,
		0.0,1.0,
	};
	double I3[] = 
	{
		1.0,0.0,0.0,
		0.0,1.0,0.0,
		0.0,0.0,1.0
	};


	M1=cria_matriz(3,3);
	copia_matriz(3,3,A,M1);
	printf("\nM1=\n");
	print_matriz(3,3,M1);

	M2=cria_transposta(3,3,M1);
	printf("\nM2=\n");
	print_matriz(3,3,M2);

	printf("\n");
	printf("\nB=\n");
	print_matriz(3,2,B);
	M3=cria_transposta(3,2,B);

	printf("\nM3=\n");
	print_matriz(2,3,M3);

	M4=cria_mult(2,2,2,E,I2);
	printf("sd\n");
	print_matriz(2,2,M4);
	printf("\n");
	M5=cria_mult(3,3,3,A,I3);
	printf("sd\n");
	print_matriz(3,3,M5);
	//print_vetor(M2,3*3);
	return 0;
}
