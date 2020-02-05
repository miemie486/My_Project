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
    start=clock();
    int n,i,nn=pow(10,8);
    double h,k1,k2,yplus,y,error;
    
    for(n=10;n<=nn;n=n*10)
    {
        h=1.0/n;
        y=1.0;
        for(i=0;i<=n;i++)
        {
            k1=df(i*h,y);
            k2=df(i*h+(3/4.0)*h,y+(3/4.0)*k1*h);
            yplus=y+((1/3.0)*k1+(2/3.0)*k2)*h;
            y=yplus;
        }
        error=fabs((f(1+h)-y)/f(1+h));
        printf("%.0f %.16f\n",log10(n),error);
    }
    end=clock();
    printf("%f\n",(end-start)/(CLOCKS_PER_SEC));
    
}