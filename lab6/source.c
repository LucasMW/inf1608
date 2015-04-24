#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"matrixSolving.h"
//#define DEBUG //uncomment this line if you want to see the computing proccess
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
			
			A[MIndex(i,j,n)]=val;
		}
	}
	A[MIndex(0,0,n)]=1;
	A[MIndex(0,1,n)]=0;
	A[MIndex(n-1,n-1,n)]=1;
	A[MIndex(n-1,n-2,n)]=0;
	A[MIndex(n-2,n-1,n)]=delta(x,n-2);
	#ifdef DEBUG
	MPrint(A,n);
	#endif

	for(i=1;i<n-1;i++) //fill vector
		b[i]=3*(delta(y,i)/delta(x,i)-delta(y,i-1)/delta(x,i-1));
	b[0]=0;
	b[n-1]=0;
	#ifdef DEBUG
	VPrint(b,n);
	#endif
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
	#ifdef DEBUG
	printf("x: %lf ",x);
	#endif
	for(i=0;i<s->n-1;i++)
	{
		//printf("%lf %lf\n ",s->x[i], s->x[i+1]);
		if(x>=s->x[i]&&x<=s->x[i+1])
		{
			#ifdef DEBUG
			printf("matched %d ",i);
			#endif
			break;
		}

	}
	#ifdef DEBUG
	printf("breaked in %d\n",i);
	#endif
	pot =x-s->x[i];
	return s->a[i]+s->b[i]*pot+ s->c[i]*(pot*pot)+s->d[i]*(pot*pot*pot);
	
}


int main (void)
{
	double x[]={1,2,3,4,5,6,7,8};
	double y[]={1,6,3,4,2,4,7,2};
	double* nx;
	double* ny;
	int n=8;
	int i;
	double d;
	Spline* s;
	FILE* output;
	s=SplineCreate(n,x,y);
	for(i=0;i<n-1;i++)
	{
		//printf("%d <%lf;%lf;%lf;%lf>\n",i,s->a[i],s->b[i],s->c[i],s->d[i]);
	}
	output=fopen("splineEval.txt","wt");
	//print values to file

	
	for(i=11;i<n*10;i++)
	{
		d=i/10.0;
		printf("evals: %16g for x %lf\n",SplineEval(s,d),d);

		fprintf(output,"%lf %lf\n",d,SplineEval(s,d));
	}
	printf("See logFile <SplineEval.txt> for exact coordinates\n");
	
	return 0;
}
