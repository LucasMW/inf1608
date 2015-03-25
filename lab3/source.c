#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
Newton-Rapson
x0 estimativa inicial
f função cuja raiz se deseja encontrar
fl derivada de f
p precisão da raiz
r endereço da variável que receberá a raiz por referência
*/
int NewtonRaphson(double x0, double(*f)(double x),double (*fl) (double x), int p, double* r)
{
	int iterations=0;
	double x,fx,flx;
	double error;
	x=x0;
	error = 0.5*pow(10.0,-1.0*p);
	fx=f(x);
	flx=fl(x);
	while(!(fabs(fx)<error))
	{
		if(fabs(flx)<0.1) //too near zero
			return 0; // Failed
		fx=f(x);
		flx=fl(x);
		x=x-fx/flx;
		printf("x:%lf\ti:%d\n",x,iterations);
		iterations++;
	}
	//Means Success
	*r=x; 
	return iterations-1;


}
double function1(double x)
{
	double sen3;
	double pot,x2;
	double result;
	sen3=sin(x);
	sen3*=sen3*sen3; //sin cubed
	x2=x*x; //x squared
	pot=x2*x; //x cubed
	result=-1;
	result-=pot;
	pot*=x;
	result-=2*pot;
	pot*=x2;
	result+=pot;
	result+=exp(sen3);
	return result;
}
double function1Derivative(double x)
{
	double sin1,sin2;
	double pot,x2;
	double result;
	x2=x*x;
	pot=x2;
	result=-3.0*pot;
	pot*=x;
	result-=8*pot;
	pot*=x2;
	result+=6*pot;
	sin1=sin(x);
	sin2=sin1*sin1; //sin squared
	result+=3*sin2*cos(x)*exp(sin2*sin1); 
	return result;
}
int IQI(double x0,double x1,double (*f)(double x), int p, double *r)
{
	int iterations=0;
	double x,xAnt,xAntAnt,xNew;
	double fx,fxAnt,fxAntAnt;
	double error;
	xAntAnt=x0;
	xAnt=x1;
	x=(x1+x0)/2.0;

	error = 0.5*pow(10.0,-1.0*p);
	fx=f(x);
	fxAnt=f(xAnt);
	fxAntAnt=f(xAntAnt);
	while(!(fabs(fx)<error))
	{

		printf("x:%16g,%d\n",x,iterations);
		sleep(1);
		//Salvar os valores das funçoes
		fx=f(x);
		fxAnt=f(xAnt);
		fxAntAnt=f(xAntAnt);
		//Recorrencia
		xNew=fxAnt*fx/(fxAntAnt-fxAnt)*(fxAntAnt-fx)*xAntAnt;
		xNew+=fxAntAnt*fx/(fxAnt-fxAntAnt)*(fxAnt-fx)*xAnt;
		xNew+=fxAntAnt*fxAnt/(fx-fxAntAnt)*(fx-fxAnt)*x;
		//Prepara as variáveis para a próxima iteração
		xAntAnt=xAnt;
		xAnt=x;
		x=xNew;

		iterations++;
	}
	//Means Success
	*r=x; 
	return iterations-1;
}
double function2(double x)
{
	return exp(x)+sin(x)-4.0;
}
int main (void)
{
	double r1,r2,r3;
	int it1,it2,it3;
	double resp;
	printf("part1\n");
	it1=NewtonRaphson(-2.0,function1,function1Derivative,6,&r1);
	it2=NewtonRaphson(1.0,function1,function1Derivative,6,&r2);
	it3=NewtonRaphson(2.0,function1,function1Derivative,6,&r3);

	printf("%16g %16g %16g\n",r1,r2,r3);
	printf("it1 %d \t it2 %d \t it3 %d\n",it1,it2,it3);

	printf("part2\n");
	IQI(0.5,2.5,function2,6,&resp)>0?printf("Answer: %16g\n",resp):printf("A resposta não pode ser obtida\n");

	return 0;
}
