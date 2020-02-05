#include<stdio.h>
#include<math.h>
#include "lapacke.h"
#define N 400
#define max -17   //动量的上限
//动量的取值范围大概在10^-19——10^-17
double h=35;
double M0=1;
double c=3*pow(10,8);
double C=30;
double g=10;

double f_(double q)
{
    double a;
    a=exp(-pow(q,4)/pow(h,4));
    return a;
}
double v(double q_out,double q_in)
{
    double av;
    av=g*f_(q_in)*f_(q_out);
    return av;
}
double f(double k_,double q)
{
    double a;
    if(k_==q)
    {
        a=M0*exp(-2*pow(k_/h,4));
    }
    else
    {
        a=(pow(k_,2)*M0*exp(-2*pow(k_/h,4))-pow(q,2)*M0*exp(-2*pow(q/h,4)))/(pow(q,2)-pow(k_,2));
    }
    
    return a;
}
int main()
{
    int i,j,kk;
    double q,q_,a1,a2,c1,c2,h_,d,hh;



            q=2.0;     //这是我们要算的R(q，q)的q值
            q_=q;
            d=0;
            for(j=1;j<=N;j++)
            {
                h_=2.0/(N);
                hh=C*M_PI*h_/(4*pow(cos(M_PI*(j*h_)/4),2));
                d=d+f(C*tan(M_PI*(j*h_)/4),q)*hh;
            }
            printf("%.10e\n",v(q_,q)/(1-g*d));
            //printf("%e\n",d);
        

    }
    
