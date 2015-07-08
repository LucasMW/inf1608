#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif
//#define DEBUG //uncomment this line if you want to see the computing proccess

#define G 9.81
#define L 10


double VerletIntegration(double t0,double t1,double y0,double v0,double (*f) (double t, double y)) //h = timestep
{
	double yAnt, y, yNext;
	double h, t;
	h=(t1-t0)/10;
	//initialization
	yAnt = y0;
	y= yAnt + h*v0;


	for(t=t0;t<t1;t+=h)
	{
		yNext= y + (y - yAnt) + f(t,y) * h *h;

		yAnt= y;
		y = yNext;
	}
	return y;
}

double calculateW(double t,double y,double h,double w,double (*f) (double t, double y))
{
	double k1,k2,k3,k4;
		k1=h*f(t,y);
		k2=h*f(t+h/2,y+k1/2.0);
		k3=h*f(t+h/2,y+k2/2.0);
		k4=h*f(t+h,y+k3);
		return w + (k1+2*(k2+k3)+k4)/6.0;

}
double AdaptativeRungeKutta(double t0, double y0, double w0, double h0, double t1, double (*f) (double t, double y), double emax)
{
	double y,y1,y2;
	double s1,s2,s3,s4;
	double k1,k2,k3,k4;
	double w,w1,w2;
	double t;
	double fty;
	double h,hNew;
	double error;
	//t1=t1+h/2.0; //interval security
	for(t=t0,y=y0,h=h0;t<t1;)
	{
		//control eval
		k1=h*f(t,y);
		k2=h*f(t+h/2,y+k1/2.0);
		k3=h*f(t+h/2,y+k2/2.0);
		k4=h*f(t+h,y+k3);
		w1 = w + (k1+2*(k2+k3)+k4)/6.0;


		s1=h*f(t,w1);
		s2=h*f(t+h/2,w+s1/2.0);
		s3=h*f(t+h/2,w+s2/2.0);
		s4=h*f(t+h,w+s3);
		y1 = y + (s1+2*(s2+s3)+s4)/6.0;



		//test eval
		s1=h*f(t,y);
		s2=h*f(t+h/4,y+s1/2.0);
		s3=h*f(t+h/4,y+s2/2.0);
		s4=h*f(t+h/2,y+s3);
		y2 = y + (s1+2*(s2+s3)+s4)/6.0;

		s1=h*f(t,y2);
		s2=h*f(t+h/4,y2+s1/2.0);
		s3=h*f(t+h/4,y2+s2/2.0);
		s4=h*f(t+h/2,y2+s3);
		y2 = y2 + (s1+2*(s2+s3)+s4)/6.0;
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

			hNew=h*pow(emax/error,-5); //5th root  
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
double RungeKutta(double t0, double t1, double h, double y0, double (*f) (double t, double y))
{
	double s1,s2,s3,s4;
	double y;
	double delta;
	double t;
	for(t=t0,y=y0;t<t1;t+=h)
	{
		if(t+h>t1)
			h=h/3.0;
		s1=h*f(t,y);
		s2=h*f(t+h/2,y+s1/2.0);
		s3=h*f(t+h/2,y+s2/2.0);
		s4=h*f(t+h,y+s3);
		y = y + (s1+2*(s2+s3)+s4)/6.0;
	}
	return y;
	
}



double FunctionSecondDerivative(double t, double theta)
{
	return -G/L * sin(theta);
}
double FunctionDerivative(double theta)
{
	return  G/L * cos(theta);
}
double Function(double t, double theta)
{
	return -G/L * sin(theta);
}
double FunctionRight(double t,double theta0)
{
	return theta0 * cos(sqrt(G/L * t));
}
int main (void)
{
	double theta;
	double theta0= M_PI/30;
	double T = 2 * M_PI * sqrt(L /G);
	double t;
	double w;
	for(t=0;t<T;t += 0.01)
	{
		printf("theta(%lf): %lf VS %lf\n",t,FunctionRight(t,theta0),VerletIntegration(0,t,theta0,0,FunctionSecondDerivative));

	}
	
	
	return 0;

}