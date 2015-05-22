#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DEBUG //uncomment this line if you want to see tha computing proccess
double Simpson (double a, double b, double (*f) (double x))
{
	double mid=(a+b)/2.0;
	double h = b-a;
	double result;
	double error;
	result=f(a)+4*f(mid)+f(b);
	result*=h/6.0;
	error= h/2.0;
	error*=error;
	error*=error;
	error/=-90.0;
	#ifdef DEBUG
	printf("error is %lf\n",error);
	#endif
	return result;
}

double DoubleSimpson(double a, double b, double (*f) (double x))
{
	double mid=(a+b)/2.0;
	double h = b-a;
	double result;
	double error;
	result=Simpson(a,mid,f)+Simpson(mid,b,f);
	error=15.0*(Simpson(a,b,f)-Simpson(a,mid,f)-Simpson(mid,b,f));
	#ifdef DEBUG
	printf("error is %lf\n",error);
	#endif
	return result;
}
double Function1(double x)
{
	return x/sqrt(x*x*x*x+1);
}
double Function2(double x)
{
	return x*x*sin(x);
}
double AdaptativeSimpson(double a, double b, double (*f) (double x), double tol)
{
	double S;
	double E;
	double D;
	double h;

	h = b-a;

	S=Simpson(a,b,f);
	D=Simpson(a,a+h/2,f)+Simpson(a+h/2,b,f);
	E = (S-D)/15.0;
	#ifdef DEBUG 
	printf("a:%lf b:%lf tol:%lf\n",a,b,tol);
	#endif	
	if(fabs(E)>tol)
	{
		#ifdef DEBUG
		printf("Entered Tol S %lf E %lf D %lf tol %lf E-tol %lf\n",S,E,D,tol,E-tol);
		#endif
		return AdaptativeSimpson(a,a+h/2.0,f,tol/2.0)+AdaptativeSimpson(a+h/2.0,b,f,tol/2.0);
	}
	else
	{
		#ifdef DEBUG
		printf("Passed S %lf E %lf D %lf tol %lf\n",S,E,D,tol);	
		#endif
		return S-E;
	}

}
static double SimpsonPlusEvalf(double a, double b, double S, double fa, double fb, double fc, int bottom, double (*f)(double x),double tol)
{
double c = (a + b)/2, h = b - a;                                                                  
  double d = (a + c)/2, e = (c + b)/2;                                                              
  double fd = f(d), fe = f(e);                                                                      
  double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
  double Sright = (h/12)*(fc + 4*fe + fb);                                                          
  double S2 = Sleft + Sright;    
  #ifdef DEBUG
  printf("jdakdn\n");
  #endif                                                                   
  if (fabs(S2 - S) <= 15*tol)                                                    
    return S2 + (S2 - S)/15;                                                                        
  return SimpsonPlusEvalf( a, c,  Sleft,  fa, fc, fd, bottom-1, f, tol/2.0) + SimpsonPlusEvalf( c, b, Sright, fc, fb, fe, bottom-1, f , tol/2.0); 
}
double AdaptativeSimpsonPlus(double a, double b, double (*f)(double x), double tol)
{
	double Sab; 
	double Sac; 
	double Scb;
	double fa;
	double fb;
	double fc;
	double c=(b-a)/2.0;
	fa=f(a);
	fb=f(b);
	fc=f(c);
	Sab=Simpson(a,b,f);
	Sac=Simpson(a,c,f);
	Scb=Simpson(c,b,f);

	return SimpsonPlusEvalf(a,b,Sac,fa,fb,fc,1000,f,tol);
}


int main (void)
{
	printf("Simpson Tests\n");
	printf("int(f1,%lf,%lf)=%lf\n",0.0,1.0,Simpson(0.0,1.0,Function1));
	printf("int(f2,%lf,%lf)=%lf\n",0.0,M_PI,Simpson(0.0,M_PI,Function2));
	printf("Double Simpson Tests\n");
	printf("int(f1,%lf,%lf)=%lf\n",0.0,1.0,DoubleSimpson(0.0,1.0,Function1));
	printf("int(f2,%lf,%lf)=%lf\n",0.0,M_PI,DoubleSimpson(0.0,M_PI,Function2));
	printf("Adaptative Simpson Tests\n");
	printf("int(f1,%lf,%lf)=%lf\n",0.0,1.0,AdaptativeSimpson(0.0,1.0,Function1,0.5*pow(10.0,-5.0)));
	printf("int(f2,%lf,%lf)=%lf\n",0.0,M_PI,AdaptativeSimpson(0.0,M_PI,Function2,0.5*pow(10.0,-5.0)));
	printf("Adaptative Simpson Plus Tests\n");
	printf("int(f1,%lf,%lf)=%lf\n",0.0,1.0,AdaptativeSimpsonPlus(0.0,1.0,Function1,0.5*pow(10.0,-5.0)));
	printf("int(f2,%lf,%lf)=%lf\n",0.0,M_PI,AdaptativeSimpsonPlus(0.0,M_PI,Function2,0.5*pow(10.0,-5.0)));

	return 0;
}
