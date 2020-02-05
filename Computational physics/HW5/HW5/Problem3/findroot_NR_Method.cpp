#include<stdio.h>
#include<math.h>

double f(double x)
{
    double a;
    a=924*pow(x,6)-2772*pow(x,5)+3150*pow(x,4)-1680*pow(x,3)+420*x*x-42*x+1;
    return a;
}
double f1(double x)
{
    double a;
    a=924*6*pow(x,5)-2772*5*pow(x,4)+3150*4*pow(x,3)-1680*3*pow(x,2)+420*2*x-42;
    return a;
}
int main()
{
    double n=1000,xn,xplus,xnn,i;
    for(i=0;i<=n;i++)
    {
        xn=i/n;
        xplus=(i+1)/n;
        xnn=xn;
        if(f(xn)*f(xplus)<0)
        {
            while(abs(xnn-xplus)>=pow(10,-13))
            {
                xnn=xn;
                xplus=xn-f(xn)/f1(xn);
                xn=xplus;
            }
            printf("%.13f\n",xn);
        }

    }
}