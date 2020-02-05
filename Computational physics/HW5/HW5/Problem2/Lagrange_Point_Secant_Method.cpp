#include<stdio.h>
#include<math.h>

double G=6.674*pow(10,-11),M=5.974*pow(10,24),m=7.348*pow(10,22),R=3.844*pow(10,8),omega=2.662*pow(10,-6);

double f(double r)
{
    double a;
    a=G*M/(r*r)-G*m/((R-r)*(R-r))-omega*omega*r;
    return a;
}

int main()
{
    double xn_1=1,xn=2,xnn=2,xplus=3;
    while(abs(xnn-xplus)>=pow(10,-11))
    {
        xnn=xn;
        xplus=xn-f(xn)*(xn-xn_1)/(f(xn)-f(xn_1));
        xn_1=xn;
        xn=xplus;
    }
    printf("%.10f",xplus);
}