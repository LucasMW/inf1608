#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define M_PI 3.14159265358979323846
double myCos(double x)
{//Taylor grau5
    double pot;
    double result;
    
    result=1;
    pot=x*x;
    result-=pot/2.0;
    pot*=pot;
    result+=pot/24.0;
    return result;

}
double residuo(double x)
{
    return -1.0/720.0*x*x*x*x*x*x;
}
void myCosTest(void)
{
    int i;
    double r,x,res;
    srand(time(NULL)); //changes random seed
    for(i=0;i<1000;i++)
    {	
        r=(double)rand()/RAND_MAX;
	//if(i%100==0)
	//	printf("%lf\n",r); //Test random purposes
        x=-M_PI/4.0+ r*M_PI/2.0;
        if(!fabs(myCos(x)-cos(x))<residuo(x))
        {
            printf("Error!\n");
            break;
        }
    }
    printf("All right\n");
}
void doubleTest(void)
{
	printf("1.2-1.0-0.2=%16g\n",1.2-1.0-0.2);
}

int main(void)
{
    printf("INF1608 ex1\n");
    printf("part1\n");
    doubleTest();
    printf("part2\n");

    
    myCosTest();
    return 0;
} 
