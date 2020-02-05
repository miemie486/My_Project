#include<stdio.h>
#include<math.h>
#include<time.h>
double df(double x,double y)
{
    double a;
    a=1-x+4*y;
    return a;
}
double f(double x)
{
    double b;
    b=x/4.0-3/16.0+(19/16.0)*exp(4*x);
    return b;
}
int main()
{
    double start,end;
    int n,i,nn=pow(10,8);
    double h,k1,k2,k3,yplus,y,error;
    double a1=1/6.0,a2=4/6.0,a3=1/6.0;
    start=clock();
    for(n=10;n<=nn;n=n*10)
    {
        h=1.0/n;
        y=1.0;
        for(i=0;i<=n;i++)
        {
            k1=df(i*h,y);
            k2=df(i*h+0.5*h,y+0.5*k1*h);
            yplus=y+k2*h;
            y=yplus;
        }
        error=fabs((f(1+h)-y)/f(1+h));
        printf("%.0f %.16f\n",log10(n),error);
    }
    end=clock();
    //printf("%f\n",(end-start)/(CLOCKS_PER_SEC));
    
}