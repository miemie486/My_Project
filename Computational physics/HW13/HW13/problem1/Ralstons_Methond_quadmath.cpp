#include<stdio.h>
#include<math.h>
#include<quadmath.h>
#include<time.h>
__float128 df(__float128 x,__float128 y)
{
    __float128 a;
    a=1-x+4.0q*y;
    return a;
}
__float128 f(__float128 x)
{
    __float128 b;
    b=x/4.0q-3/16.0q+(19/16.0q)*expq(4.0q*x);
    return b;
}
int main()
{
    double start,end;
    start=clock();
    int n,i,nn=pow(10,8);
    __float128 h,k1,k2,yplus,y,error;
    for(n=10;n<=nn;n=n*10)
    {
        h=1.0q/n;
        y=1.0q;
        for(i=0;i<=n;i++)
        {
            k1=df(i*h,y);
            k2=df(i*h+(3/4.0q)*h,y+(3/4.0q)*k1*h);
            yplus=y+((1/3.0q)*k1+(2/3.0q)*k2)*h;
            y=yplus;
        }
        error=fabsq((f(1+h)-y)/f(1+h));
        printf("%.0Qf %.32Qf\n",log10q(n),error);
    }
    end=clock();
    printf("%f\n",(end-start)/(CLOCKS_PER_SEC));
    
}