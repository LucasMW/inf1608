#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define DEBUG //uncomment this line if you want to see the computing proccess

#ifndef Sleep
#include <windows.h>
#endif


double AdaptativeEuler(double t0, double y0, double h0, double t1, double (*f) (double t, double y), double emax)
{
	double y,y1,y2;
	double t;
	double fty;
	double h,hNew;
	double error;
	//t1=t1+h/2.0; //interval security
	for(t=t0,y=y0,h=h0;t<t1;)
	{
		y1=y+h*f(t,y);
		y2=y+h/2.0*f(t,y);
		y2=y2+h/2.0*f(t,y2); 
		error=fabs(y2-y1);
		
		if(t+h>t1)
		{
				printf("t %lf will pass t1 %lf with h %lf\n",t,t1,h);
				h=t1-t;
				printf("new h %lf\n",h);
		}
		else
		{	if(error<emax)
			{
				printf("step h %lf accepted y: %lf",h,y);
				t+=h;
				y=y2;
				printf(" becoming y %lf and t %lf\n",y,t);
			}
			else
			{
				printf("error %lf greater tham emax %lf\n",error,emax);
			}
			hNew=h*sqrt(emax/error);
			if(hNew>1.2*h)
			hNew=1.2*h;
			h=hNew;
			printf("h is now %lf\n",h);
#ifdef DEBUG
			Sleep(1000);
#endif
		}
		
		
	}
	return y;
}
double FunctionDerivative(double t, double y)
{
	return (y+t*t)*t;
}
double FunctionRight(double t)
{
	return exp(t*t/2.0)-t*t-2;
}
int main (void)
{
	double error;
	double r,y;
	double t=2.4;
	double t0=0;
	double y0=-1.0;
	double h0=0.001;
	double emax= 0.0001;
	r = AdaptativeEuler(t0,y0,h0,t,FunctionDerivative,emax);
	y=FunctionRight(t);
	error=(r-y)/y;

	printf("r %lf vs y %lf relative error : %lf \n",r,y,error);
		

	return 0;

}
