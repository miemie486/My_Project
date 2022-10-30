#include<stdio.h>
#include<math.h>

double f(double x)
{
    double a;
    a=924*pow(x,6)-2772*pow(x,5)+3150*pow(x,4)-1680*pow(x,3)+420*x*x-42*x+1;
    return a;
}

int main()
{
    double n=1000,xn_1,xn,xplus,xnn,i;
    for(i=0;i<=n;i++)
    {
        xn=(i+1)/n;
        xn_1=i/n;
        xplus=(i+2)/n;
        xnn=xn;
        if(f(xn)*f(xplus)<0)
        {
            while(abs(xnn-xplus)>=pow(10,-13))
            {
                xnn=xn;
                xplus=xn-f(xn)*(xn-xn_1)/(f(xn)-f(xn_1));
                xn_1=xn;
                xn=xplus;
            }
            printf("%.13f\n",xn);
        }

    }
}