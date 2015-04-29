#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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
			T[index(i,j,n)]=A[index(j,i,n)];
	return T;
}
double* cria_mult(int m,int n,int l, double* A, double* B)
{
	double * C;
	int i,j,k;
	C=cria_matriz(m,n);
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			for(k=0;k<n;k++)
				C[index(i,j,n)]=A[index(i,k,n)]*B[index(k,j,n)];
}

double* cria_mmq(int m,int n,double* A,double* b)
{
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
	M1=cria_matriz(3,3);
	copia_matriz(3,3,A,M1);

	print_matriz(3,3,M1);

	M2=cria_transposta(3,3,M1);
	print_matriz(3,3,M2);

	printf("\n");
	print_matriz(3,2,B);
	M3=cria_transposta(3,2,B);
	print_matriz(2,3,M3);
	//print_vetor(M2,3*3);
	return 0;
}
