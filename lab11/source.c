#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>
//#define DEBUG //uncomment this line if you want to see the computing proccess

#ifndef M_PI
#define M_PI 3.14159265358979323846
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
		//control eval
		y1=y+h*f(t,y);

		//test eval
		y2=y+h/2.0*f(t,y);
		y2=y2+h/2.0*f(t+h/2.0,y2); 

		//step error
		error=fabs(y2-y1);
		
		if(t+h>t1)
		{
				#ifdef DEBUG
				printf("t %lf will pass t1 %lf with h %lf\n",t,t1,h);
				#endif
				h=t1-t;
				#ifdef DEBUG
				printf("new h %lf\n",h);
				#endif
		}
		else
		{	
			if(error<emax) //error is tolerable
			{
				#ifdef DEBUG
				printf("step h %lf accepted y: %lf",h,y);
				#endif
				t+=h; //accept step
				y=y2; //accept step value
				#ifdef DEBUG
				printf(" becoming y %lf and t %lf\n",y,t);
				#endif
			}
			else
			{
				#ifdef DEBUG
				printf("error %lf greater tham emax %lf\n",error,emax);
				#endif
			}

			hNew=h*sqrt(emax/error);
			if(hNew > 1.2*h)
			{
				hNew = 1.2*h;
			}
			h=hNew;
			#ifdef DEBUG
			printf("h is now %lf\n",h);
			#endif

		}
		
		
	}
	return y;
}
double TimeToY1(double t0, double y0, double y1, double h, double (*f) (double t, double y))
{

	double s1,s2,s3,s4;
	double y;
	double delta;
	double t;
	//t1=t1+h/2.0; //interval security
	for(t=t0,y=y0; y<y1;t+=h)
	{
		// if(t+h>t1)
		// 	h=h/3.0;
		s1=h*f(t,y);
		s2=h*f(t+h/2,y+s1/2.0);
		s3=h*f(t+h/2,y+s2/2.0);
		s4=h*f(t+h,y+s3);
		y = y + (s1+2*(s2+s3)+s4)/6.0;
		
		#ifdef DEBUG
		printf("RungeKutta t: %lf y:%lf \n",t,y);
		#endif
	}

	return t;
	

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
	double error,rError;
	double r,y;
	double t=2.4;
	double t0=0;
	double y0=-1.0;
	double h0=0.001;
	double emax= 0.0001;
	r = AdaptativeEuler(t0,y0,h0,t,FunctionDerivative,emax);
	y=FunctionRight(t);
	error=fabs(r-y);
	rError=error/y;

	printf("r %lf vs y %lf relative error : %lf  absolute error %lf\n",r,y,rError,error);
	r=TimeToY1(t0,y0,y,h0,FunctionDerivative);
	printf("r %lf vs t %lf absolute error %lf\n",r,t,fabs(r-t));


	return 0;

}
