#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"matrixSolving.h"

struct spline //Cubic Splines
{	
	int n; //number of points
	double* x; //vector of x coordinates
	double* y; //vector of y coordinates
	double* a; //vector of a coeficients
	double* b; //vector of b coeficients
	double* c; //vector of c coeficients
	double* d; //vector of d coeficients
};

typedef struct spline Spline;

void VectorCopy(double* org, double* dest,int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		dest[i]=org[i];
	}

}
static double delta(double* v,int i)
{
	return v[i+1]-v[i];
}
Spline* SplineCreate (int n,double* x, double* y)
{
	int i,j;
	double val;
	Matrix A;
	double* b;
	Spline* spline;
	//memory allocating
	A=(Matrix)malloc(sizeof(double)*n*n);
	b=(double*)malloc(sizeof(double)*n);
	spline=(Spline*)malloc(sizeof(Spline));
	spline->a=(double*)malloc(sizeof(double)*n);
	spline->b=(double*)malloc(sizeof(double)*n);
	spline->c=(double*)malloc(sizeof(double)*n);
	spline->d=(double*)malloc(sizeof(double)*n);
	spline->x=(double*)malloc(sizeof(double)*n);
	spline->y=(double*)malloc(sizeof(double)*n);
	spline->n=n;
	VectorCopy(x,spline->x,n);
	VectorCopy(y,spline->y,n);
	for(i=0;i<n;i++) //fill matrix
	{
		for(j=0;j<n;j++)
		{
			val=0;
			if(i>1&&j>1&&i<n-1,j<n-1)
			{
				if(i==j) //diagonal
				{
					val=2*(delta(x,i-1)+delta(x,i));
				}
				else if(j==(i+1)) //diagonal direita
				{
					val=delta(x,i);
				}
				else if(j==(i-1)) //diagonal esquerda
				{
					val=delta(x,i-1);
				}
			}
			A[MIndex(i,j,n)]=val;
		}
	}
	A[MIndex(0,0,n)]=1;
	A[MIndex(0,1,n)]=0;
	A[MIndex(n-1,n-1,n)]=1;
	A[MIndex(n-1,n-2,n)]=0;
	A[MIndex(n-2,n-1,n)]=delta(x,n-2);
	MPrint(A,n);

	for(i=1;i<n-1;i++) //fill vector
		b[i]=3*(delta(y,i)/delta(x,i)-delta(y,i-1)/delta(x,i-1));
	b[0]=0;
	b[n-1]=0;
	VPrint(b,n);
	NaiveGaussElimination(n,A,b);
	BackSubstitution(n,A,b,spline->c); //solve system
	for(i=0;i<n-1;i++) //fill spline coeficients
	{
		spline->a[i]=y[i];
		spline->b[i]=delta(y,i)/delta(x,i)-delta(x,i)/3.0*(2*spline->c[i]+spline->c[i+1]);
		spline->d[i]=(spline->c[i+1]-spline->c[i])/3*delta(x,i);
	}
	return spline;
}
double SplineEval(Spline* s,double x)
{
	int i;
	double pot;
	for(i=0;i<s->n-1;i++)
	{
		if(x>s->x[i]&&x<s->x[i+1])
		{
			printf("matched %d ",i);
			break;
		}
	}
	printf("breaked in %d\n",i);
	pot =x-s->x[i];
	return s->a[i]+s->b[i]*pot+ s->c[i]*(pot*pot)+s->d[i]*(pot*pot*pot);
	
}
int main (void)
{
	double x[]={1,2,3,4,5,6,7,8};
	double y[]={1,6,3,4,2,4,7,2};
	int n=8;
	int i;
	double d;
	Spline* s;
	s=SplineCreate(n,x,y);
	for(i=0;i<n;i++)
	{
		printf("%d <%lf;%lf>\n",i,s->c[i],s->d[i]);
	}
	for(i=11;i<n*10;i++)
	{
		d=i/10.0;
		printf("i: %d evals: %16g d:%lf\n",i,SplineEval(s,d),d);
	}
	
	return 0;
}
